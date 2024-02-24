using LinearAlgebra
using Dates
using Printf
using Statistics
using FITSIO
using Snap
using Glob
using Interpolations: Interpolations
using CoordinateTransformations: CoordinateTransformations
using StaticArrays
using OffsetArrays: OffsetArrays
using ImageTransformations: ImageTransformations
using ProgressLogging
using AstroImages

"""
This function peforms sky background calibration, with sky frames
and/or chopping.

It must be used on calibrated images that may already have a coarse BG subtraction.
It finds the star and then masks it, and then re-performs calibration from raws
including a better BG sub (due to the masking) and supports chopping BG sub (auto-
detected if star has moved a lot during the sequence).
"""
function bgsub(conf_fname; verbose=true, savedata=true, showplots=true, force=false)
    println(Snap.banner()*"\n\nBGSUB")
   
    # Read the configuration
    conf = Snap.readconfig(conf_fname)
    # cal = conf["calibrate"]
    # λoverD = conf["telescope"]["λoverD"];

    profile = replace(conf_fname, ".toml" => "")
    profilename = splitpath(profile)[end-1]
    try
        global profilename = replace(joinpath(splitpath(profile)[end-1:end]...), r"\W"=>"-")
    catch
    end
    profdir = dirname(profile)
    cd(profdir)
    profdir= ""
    caldir = joinpath(profdir,"cal",basename(profile))
    if !isdir(joinpath(profdir,"cal"))
        savedata && mkpath(joinpath(profdir,"cal"))
    end

    sky_psf_radius_px = get(conf["calibrate"], "sky_psf_radius_px", 75.0)

    @info "Loading skys"
    # dark = load("$caldir.dark.fits.gz")
    # darkpsf = load("$caldir.darkpsf.fits.gz")
    # skypsf = load("$caldir.skypsf.fits.gz")
    # flat = load("$caldir.flat.fits.gz")
    skys = load("$caldir.sky.fits.gz", :)
    # bpm = load("$caldir.bpm.fits.gz")

    fnames_0 = Glob.glob("$caldir.*.cal.fits.gz")
    @info "Loading lights"
    cals = AstroImageMat[]
    outfnames = String[]
    for fname in fnames_0
        outfname = replace(fname, ".fits"=>".bgsub.fits",)
        # Load all images if chopping. Otherwise only load if changed (and not forced)
        if get(conf["calibrate"], "chopped", true) || 
            force ||
            (isfile(outfname) && Base.Filesystem.mtime(outfname) > Base.Filesystem.mtime(fnames[i]))
            push!(cals, load(fname))
            push!(outfnames, outfname)
            continue
        end
    end
    if isempty(cals)
        return
    end

    @info "Masking star"
    # Note, this is multi-threaded so afterwards cals_masked are out of order.
    cals_masked = Any[nothing for _ in cals]
    Threads.@threads :dynamic for i in eachindex(cals)
        cal = cals[i]
        # checknansdebug(cal,fname)
        calmasked = deepcopy(cal)
        xs = axes(calmasked,1) .- calmasked["STAR-X"]
        ys = axes(calmasked,2) .- calmasked["STAR-Y"]
        rs = @. sqrt(xs^2 + ys'^2)
        calmasked[rs .<= sky_psf_radius_px] .= NaN
        cals_masked[i] = mapwindow(median, calmasked, (3,3))
    end

    # remove bad pixels using a median filter
    skys = map(skys) do sky
        mapwindow(median, sky[:,:], (3,3))
    end

    # We need to subract off the median values of all images for the least squares to work well.
    # Measure the median of science images from their masked counter parts (so star doesn't bias)
    meds = map(cals_masked) do cal
        median(filter(isfinite, cal))
    end
    cals = broadcast(cals, meds) do img, med
        img .- med
    end
    cals_masked = broadcast(cals_masked, meds) do img, med
        img .- med
    end
    # Do same for sky images
    meds = map(skys) do sky
        median(filter(isfinite, sky))
    end
    skys = broadcast(skys, meds) do img, med
        checknansdebug(img)
        (img .- med)[:,:]
    end


    if savedata && get(conf["calibrate"], "write_calmasked", false)
        @progress "Writing" for (i,frame) in enumerate(cals_masked)
            AstroImages.writefits(@sprintf("%s.%04d.3.cal.masked.fits.gz", caldir, i),frame)
        end
    end
    # gif_fname = "$caldir.anim.3.rawmasked4sky.gif"
    # @ANIMATE
    # frames = stack(raws_masked[begin:4:end,begin:4:end])
    # save(gif_fname, imview(frames), fps=10)


    sky_N_most_similar = get(conf["calibrate"], "sky_N_most_similar", typemax(Int))
    

    star_x_y = map(cals) do cal
        cal["STAR-X"], cal["STAR-Y"]
    end

    
    @info "Calibrating raw images (sky and/or choppping)"
    savers = []
    Threads.@threads :dynamic for i in eachindex(cals)
        outfname = outfnames[i]
        A = cals[i]

        if get(conf["calibrate"], "chopped", true)

            # Determine which frames are far enough away to use as a BG reference
            this_pos = star_x_y[i]
            distances = map(xy->sqrt((this_pos[1] - xy[1])^2 + (this_pos[2] - xy[2])^2), star_x_y)
            
            # These are the frames we will be allowed to use as bg references
            ii_allowed = findall(distances .>= sky_psf_radius_px)
            # Both regular sky frames and light frames that are far enough away to be used as sky backgrounds (chopping)
            chop_skys = vcat(cals[ii_allowed], skys)
            masked_chop_skys = vcat(cals_masked[ii_allowed], skys)
        else
            # No chopping, just do least squares fit to BG
            chop_skys = skys
            masked_chop_skys = skys
        end
        # Rank by correlation
        correlation_vector = map(masked_chop_skys) do B
            mask = isfinite.(A) .& isfinite.(B)
            if count(mask) == 0 || all(A[mask] .== B[mask])
                return 0.0
            end
            return abs(cor(A[mask],B[mask],))
        end
        II = sortperm(correlation_vector, rev=true)[1:min(end,sky_N_most_similar)]
        
        # Lets by rank time!
        # t0 = DateTime(cal.headers["DATE-OBS"]*" "*cal.headers["UTC"],dateformat"yyyy-mm-dd HH:MM:SS.s")
        # correlation_vector = map(skysmasked) do sky
        #     ti = DateTime(sky.headers["DATE-OBS"]*" "*sky.headers["UTC"],dateformat"yyyy-mm-dd HH:MM:SS.s")
        #     return abs(t0-ti)
        # end 
        # II = sortperm(correlation_vector)[1:min(end,sky_N_most_similar)]

        # Perform the calibration
        skycube = stack(chop_skys[II])
        skycubemasked = stack(masked_chop_skys[II])


        ## least squares subtraction of sky background

        # Tile the image into large chunks, doing the fitting one chunk at a time instead 
        # of finding a global solution
        newsky = deepcopy(A)
        fill!(newsky, 0)

        tile_count = get(conf["calibrate"], "bgsub-tile-count", 2)
        tiles_x = range(1, stop=size(skycubemasked,1), length=tile_count+1)
        tiles_x_start = floor.(Int, tiles_x[begin:end-1])
        tiles_x_stop = floor.(Int, tiles_x[begin+1:end])
        tiles_y = range(1, stop=size(skycubemasked,2), length=tile_count+1)
        tiles_y_start = floor.(Int, tiles_y[begin:end-1])
        tiles_y_stop = floor.(Int, tiles_y[begin+1:end])
        j = 0
        for (start_x,stop_x) in zip(tiles_x_start,tiles_x_stop), 
            (start_y,stop_y) in zip(tiles_y_start,tiles_y_stop)
            tilemask = (start_x:stop_x,start_y:stop_y)
            j += 1

            # @info "tile BG sub" j tilemask
                        
            # speed it up by reducing the image and cube dimensions 2X in spatial coorinates
            # This is a balancing act: the more it is downsized, the more it will focus on large
            # scale structure at the expense of pixel to pixel structure.
            # skycubemasked_restrict = skycubemasked
            # cal_masked_restrict = cals_masked[i]
            skycubemasked_restrict = AstroImages.restrict(skycubemasked[tilemask...,:],(1,2))
            cal_masked_restrict = AstroImages.restrict(cals_masked[i][tilemask...],(1,2))

            finite_mask = isfinite.(cal_masked_restrict) .& all(isfinite, skycubemasked_restrict ,dims=3)[:,:]
            if count(finite_mask) == 0
                @warn "empty tile, no bg subtraction" tilemask
                return
            end

            # Also ignore pixels outside a certain quantile range
            l,u = quantile(@view(cal_masked_restrict[finite_mask]), (0.001, 0.999))
            finite_mask .&= l .<= cal_masked_restrict .<= u
            l,u = quantile(@view(skycubemasked_restrict[finite_mask,:]), (0.001, 0.999))
            finite_mask .&= all(l .<= skycubemasked_restrict .<= u,dims=3)[:,:]
            finite_mask = BitMatrix(finite_mask[:,:])

            S = reshape(skycubemasked_restrict[finite_mask,:],count(finite_mask), :)
            b = vec(cal_masked_restrict[finite_mask])
            coeff = S \ b
            if any(!isfinite, coeff)
                @warn "non finite sky least square coefficients" 
                coeff .= 1 / length(coeff)
            end

            # Apply subtraction using un-masked versions
            newsky[tilemask...] = sum(skycube[tilemask...,:] .* reshape(coeff,1,1,:),dims=3)[:,:]


            ##############################################
            # Stuck debugging why the BG subtract only works when I run it on ONE
            ##################################################
            # if i == 1
            # if contains(outfname, "n0145")
            # Handy for debugging: uncomment this to dump the various cubes used to make the background
            # for this tile.
                # AstroImages.writefits(
                #     replace(outfname, ".fits"=>".A.$j.fits"),
                #     cal_masked_restrict .- sum(skycubemasked_restrict.* reshape(coeff,1,1,:),dims=3)[:,:]
                # )
                # AstroImages.writefits(
                #     replace(outfname, ".fits"=>".B.$j.fits"),
                #     cal_masked_restrict,
                #     A,
                # )
                # AstroImages.writefits(replace(outfname, ".fits"=>".C.$j.fits"), collect(skycubemasked_restrict))
                # AstroImages.writefits(replace(outfname, ".fits"=>".D.$j.fits"), UInt8.(finite_mask))
                # AstroImages.writefits(replace(outfname, ".fits"=>".E.$j.fits"), reshape(coeff,:,1))
            # end
        end
    
        corrected = copyheader(A, A .- newsky)

        push!(corrected, History, "$(Date(Dates.now())): LOCI background subtraction.")
        corrected["BGTILEN"] = tile_count


        # if !contains(outfname, "n0146")
        # end
        push!(savers, (outfname, corrected))
        # display(imview(corrected))

        # saver = @async AstroImages.writefits(@sprintf("%s.%04d.4.synthsky.fits.gz", caldir, i),newsky)
        # push!(savers, saver)
        println(outfname, "\t($i)")
    end

    for (fname, img) in savers
        AstroImages.writefits(fname,img)
    end


    return


end




# """
# ridge_regression(Y,X,λ,β₀=0)

# Calculate ridge regression estimate with target vector β₀.

# https://discourse.julialang.org/t/hello-i-am-stuck-with-a-problem-regarding-linear-regression-and-ridge-regression/75382/11?u=sefffal
# """
# function ridge_regression(Y,X,λ,β₀=0)
#     K = size(X,2)
#     isa(β₀,Number) && (β₀=fill(β₀,K))
#     local b
#     try
#         b = (X'X+λ*I)\(X'Y+λ*β₀)      #same as inv(X'X+λ*I)*(X'Y+λ*β₀)
#     catch err
#         @error "error in ridge ridge_regression" exception=(err, catch_backtrace())
#         return NaN .* (X \ Y)
#     end
#     return b
# end


function checknansdebug(img,fname="")
    if any(!isfinite, img[200:end-200,200:end-200])
        error("NANs $fname")
    end
end

bgsub_nirc2(args...)=error()