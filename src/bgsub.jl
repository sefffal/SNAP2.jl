using LinearAlgebra
using Dates
using Printf
using Statistics
using FITSIO
using Glob
using Interpolations: Interpolations
using CoordinateTransformations: CoordinateTransformations
using StaticArrays
using OffsetArrays: OffsetArrays
using ImageTransformations: ImageTransformations
using ProgressLogging
using AstroImages
using TOML
export bgsub

"""
This function peforms sky background calibration, with sky frames
and/or chopping.

It must be used on calibrated images that may already have a coarse BG subtraction.
It finds the star and then masks it, and then re-performs calibration from raws
including a better BG sub (due to the masking) and supports chopping BG sub (auto-
detected if star has moved a lot during the sequence).
"""
function bgsub(conf_fname, pattern=nothing, skypattern=nothing; verbose=true, savedata=true, showplots=true, force=false)
    println(SNAP.banner()*"\n\nBGSUB")
   
    # Read the configuration
    conf = TOML.parsefile(conf_fname)
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

    if isnothing(pattern) 
        pattern="$caldir.*.cal.fits.gz"
    end
    if isnothing(skypattern) 
        skypattern="$caldir.sky.fits.gz"
    end

    # How far apart do they have to be in order to be included as a ref?
    sky_psf_radius_px = get(conf["calibrate"], "sky_psf_radius_px", 75.0)
    # How far out should we mask? (often the same, but sometimes larger.)
    psf_mask_radius = get(conf["calibrate"], "psf_mask_radius", sky_psf_radius_px)

    @info "Loading skys"
    # dark = load("$caldir.dark.fits.gz")
    # darkpsf = load("$caldir.darkpsf.fits.gz")
    # skypsf = load("$caldir.skypsf.fits.gz")
    # flat = load("$caldir.flat.fits.gz")
    skys = load(skypattern, :)
    # bpm = load("$caldir.bpm.fits.gz")

    fnames_0 = Glob.glob(pattern)
    @info "Loading lights"
    cals = AstroImageMat[]
    outfnames = String[]
    for fname in fnames_0
        outfname = replace(fname, ".fits"=>".bgsub2.fits",)
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
    # Create a copy of each calibrated image where the star has been masked out.
    # Note, this is multi-threaded so afterwards cals_masked are out of order.
    cals_masked = Any[nothing for _ in cals]
    Threads.@threads :dynamic for i in eachindex(cals)
        cal = cals[i]
        calmasked = deepcopy(cal)
        xs = axes(calmasked,1) .- calmasked["STAR-X"]
        ys = axes(calmasked,2) .- calmasked["STAR-Y"]
        rs = @. sqrt(xs^2 + ys'^2)
        calmasked[rs .<= psf_mask_radius] .= NaN
        cals_masked[i] = copyheader(cal, mapwindow(median, calmasked, (3,3)))
    end

    # remove bad pixels using a median filter
    skys = map(skys) do sky
        copyheader(sky,mapwindow(median, sky[:,:], (3,3)))
    end

    # We need to subract off the median values of all images for the least squares to work well.
    # Measure the median of science images from their masked counterparts (so star doesn't bias)
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
    skys = broadcast(skys, meds) do sky, med
        copyheader(sky, (sky .- med)[:,:])
    end

    if savedata && get(conf["calibrate"], "write_calmasked", false)
        @progress "Writing" for (i,frame) in enumerate(cals_masked)
            AstroImages.writefits(@sprintf("%s.%04d.3.cal.masked.fits.gz", caldir, i),frame)
        end
    end
    sky_N_most_similar = get(conf["calibrate"], "sky_N_most_similar", typemax(Int))


    
    @info "Calibrating raw images (sky and/or choppping)"
    savers = []
    Threads.@threads :dynamic for i in eachindex(cals)
        outfname = outfnames[i]
        A = cals[i]
        
        # If NIRC2 and L or M band, then use Jayke's BG ADI technique on the skys
        # Use two sets of BG references: sky images, and sky images transformed to match rotator position
        if occursin("NIRC2", get(header(cals[i]), "CURRINST", "")) && cals[i]["WAVE"] > 3.0e-6
            @info "Preparing extra sky BG refs to account for K-mirror position"
            use_skys = [
                skys;
                [
                    prepare_nirc2_kmirror_rotated_bgrefs(A, sky)
                    for sky in skys
                    if haskey(sky, "ROTPPOSN")
                ]
            ]
        else
            use_skys = skys
        end

        # If doing chopping, use other calibrated science images as BG references if the star is 
        # far enough away
        if get(conf["calibrate"], "chopped", true)
            # If NIRC2 and L or M band, then use Jayke's BG ADI technique as well.
            # Use two sets of BG references: chopped science images,
            # and chopped science images transformed to match rotator position.
            if occursin("NIRC2", get(header(cals[i]), "CURRINST", "")) && cals[i]["WAVE"] > 3.0e-6
                @info "Preparing extra chopped BG refs to account for K-mirror position"
                use_cals = [cals; prepare_nirc2_kmirror_rotated_bgrefs.(Ref(A), cals)]
                use_cals_masked = [cals_masked; prepare_nirc2_kmirror_rotated_bgrefs.(Ref(A), cals_masked)]
            else
                # Otherwise proceed without the extra K mirror transformed references.
                use_cals = cals
                use_cals_masked = cals_masked
            end

            # Determine which frames are far enough away to use as a BG reference
            star_x_y = map(use_cals) do cal
                cal["STAR-X"], cal["STAR-Y"]
            end
            this_pos = star_x_y[i]
            distances = map(xy->sqrt((this_pos[1] - xy[1])^2 + (this_pos[2] - xy[2])^2), star_x_y)
            
            # These are the frames we will be allowed to use as bg references
            ii_allowed = findall(distances .>= sky_psf_radius_px)
            # Both regular sky frames and light frames that are far enough away to be used as sky backgrounds (chopping)
            chop_skys = vcat(use_cals[ii_allowed], use_skys)
            masked_chop_skys = vcat(use_cals_masked[ii_allowed], use_skys)
        else
            # If not doing chopping, just use sky images and not other calibrated science images
            chop_skys = use_skys
            masked_chop_skys = use_skys
        end
        # Rank by correlation. The user can allow us to only select N most correlated BG references
        correlation_vector = map(masked_chop_skys) do B
            mask = isfinite.(A) .& isfinite.(B)
            if count(mask) == 0 || all(A[mask] .== B[mask])
                return 0.0
            end
            return abs(cor(A[mask],B[mask],))
        end
        II = sortperm(correlation_vector, rev=true)[1:min(end,sky_N_most_similar)]
        
        # Alternative: rank most similar images purely by distance in time
        # t0 = DateTime(cal.headers["DATE-OBS"]*" "*cal.headers["UTC"],dateformat"yyyy-mm-dd HH:MM:SS.s")
        # correlation_vector = map(skysmasked) do sky
        #     ti = DateTime(sky.headers["DATE-OBS"]*" "*sky.headers["UTC"],dateformat"yyyy-mm-dd HH:MM:SS.s")
        #     return abs(t0-ti)
        # end 
        # II = sortperm(correlation_vector)[1:min(end,sky_N_most_similar)]

        # Now perform the BG calibration


        skycube = stack(chop_skys[II])
        skycube[.!isfinite.(skycube)] .= 0

        skycubemasked = stack(masked_chop_skys[II])
        skycubemasked[.!isfinite.(skycubemasked)] .= 0

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
        # Go through tiles
        for (start_x,stop_x) in zip(tiles_x_start,tiles_x_stop), 
            (start_y,stop_y) in zip(tiles_y_start,tiles_y_stop)
            tilemask = (start_x:stop_x,start_y:stop_y)
            j += 1

            # speed it up by reducing the image and cube dimensions 2X in spatial coorinates
            # This is a balancing act: the more it is downsized, the more it will focus on large
            # scale structure at the expense of pixel to pixel structure.
            # skycubemasked_restrict = skycubemasked
            # cal_masked_restrict = cals_masked[i]
            skycubemasked_restrict = AstroImages.restrict(skycubemasked[tilemask...,:],(1,2))
            cal_masked_restrict = AstroImages.restrict(cals_masked[i][tilemask...],(1,2))
            
            skycubemasked_restrict[.!isfinite.(skycubemasked_restrict)] .= 0
            cal_masked_restrict[.!isfinite.(cal_masked_restrict)] .= 0

            finite_mask = isfinite.(cal_masked_restrict) .& all(isfinite, skycubemasked_restrict ,dims=3)[:,:]
            if count(finite_mask) == 0
                @warn "empty tile, no bg subtraction" tilemask
                continue
            end

            # Also ignore pixels outside a certain quantile range
            l,u = quantile(@view(cal_masked_restrict[finite_mask]), (0.001, 0.999))
            finite_mask .&= l .<= cal_masked_restrict .<= u
            l,u = quantile(@view(skycubemasked_restrict[finite_mask,:]), (0.001, 0.999))
            finite_mask .&= all(l .<= skycubemasked_restrict .<= u,dims=3)[:,:]
            finite_mask = BitMatrix(finite_mask[:,:])

            # Now calculate the coefficients of the subtraction using least-squares
            S = reshape(skycubemasked_restrict[finite_mask,:],count(finite_mask), :)
            b = vec(cal_masked_restrict[finite_mask])
            coeff = S \ b
            if any(!isfinite, coeff)
                @warn "non finite sky least square coefficients" 
                coeff .= 1 / length(coeff)
            end

            # Apply subtraction using un-masked versions
            newsky[tilemask...] = sum(skycube[tilemask...,:] .* reshape(coeff,1,1,:),dims=3)[:,:]
        end
    
        corrected = copyheader(A, A .- newsky)

        push!(corrected, History, "$(Date(Dates.now())): LOCI background subtraction.")
        corrected["BGTILEN"] = tile_count


        # push!(savers, (outfname, corrected))
        push!(savers, (outfname, corrected))
        # AstroImages.writefits(outfname,corrected)
        println(outfname, "\t($i) calc")
    end

    for (fname, img) in savers
        AstroImages.writefits(fname,img)
        println(fname, "\t write")
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




    #=

TODO: If inst is NIRC2 and wavelength > ~3.0 micron, then copy all refs and create rotated versions
following Jayke's advice:


Prepare: cals, cals_masked, and skys
=#
using Rotations
"""
Given an reference image, rotate and transform it such that it tracks the background of a given target image.
"""
function prepare_nirc2_kmirror_rotated_bgrefs(targ::AstroImage, ref::AstroImage)

    # Rotate image keeping same axes
    # Transform STAR-X and STAR-Y headers too.

    center_dist = 13.75 # arcsec (from Jim Lyke)
    pixel_scale = 0.009942 # arcsec/pixel narrowcam
    pixel_radius = center_dist/pixel_scale
    rotpposn_targ = targ["ROTPPOSN"] - 0.7 # deg
    rotpposn_ref = ref["ROTPPOSN"] - 0.7 # deg

    # This works!
    # tform_targ = Translation((0, -pixel_radius)) ∘ LinearMap(RotMatrix(deg2rad(-rotpposn_targ/2))) #∘ 
    # warp(
    #     ImageFiltering.centered(targ), tform_targ, (-1750:1750, -1750:1750))#, (1:ceil(Int,pixel_radius), 1:ceil(Int,pixel_radius)))

    tform_targ = 
        CoordinateTransformations.recenter(
            LinearMap(RotMatrix(deg2rad(-rotpposn_targ/2))),
            (0, -pixel_radius)
        )
    tform_ref =  
        CoordinateTransformations.recenter(
            LinearMap(RotMatrix(deg2rad(-rotpposn_ref/2))),
            (0, -pixel_radius)
        )

    match_ref_to_targ_tform = inv(tform_targ) ∘ tform_ref
    ref_cen = ImageFiltering.centered(ref)
    matched = copyheader(ref, collect(warp(
        ref_cen,
        match_ref_to_targ_tform,
        # (-1700:1700, -1700:1700)
        axes(ref_cen)
        # (-1000:1000,-1000:1000)
    ))) 
    if haskey(ref, "STAR-X")
        newx, newy = inv(match_ref_to_targ_tform)((
            ref["STAR-X"]-size(ref,1)/2,
            ref["STAR-Y"]-size(ref,2)/2
        ))
        matched["STAR-X"] = newx + size(ref,1)/2
        matched["STAR-Y"] = newy + size(ref,2)/2
    end
    return matched 
end









"""
This function peforms sky background calibration, with sky frames
and/or chopping.

It must be used on calibrated images that may already have a coarse BG subtraction.
It finds the star and then masks it, and then re-performs calibration from raws
including a better BG sub (due to the masking) and supports chopping BG sub (auto-
detected if star has moved a lot during the sequence).
"""
function bgsub2(
    pattern=nothing,
    skypattern=nothing;
    # How far apart do they have to be in order to be included as a ref?
    sky_psf_radius_px=75.,
    # How far out should we mask? (often the same, but sometimes larger.)
    psf_mask_radius=sky_psf_radius_px,
    chopped=true,
    write_calmasked=false,
    sky_N_most_similar = typemax(Int),
    tile_count = 2,

    savedata=true,
    force=false,
)
    println(SNAP.banner()*"\n\nBGSUB")
   

    @info "Loading skys"
    # dark = load("$caldir.dark.fits.gz")
    # darkpsf = load("$caldir.darkpsf.fits.gz")
    # skypsf = load("$caldir.skypsf.fits.gz")
    # flat = load("$caldir.flat.fits.gz")
    skys = load(skypattern, :)
    # bpm = load("$caldir.bpm.fits.gz")

    fnames_0 = Glob.glob(pattern)
    @info "Loading lights"
    cals = AstroImageMat[]
    outfnames = String[]
    for fname in fnames_0
        outfname = replace(fname, ".fits"=>".bgsub2.fits",)
        # Load all images if chopping. Otherwise only load if changed (and not forced)
        if chopped || 
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
    # Create a copy of each calibrated image where the star has been masked out.
    # Note, this is multi-threaded so afterwards cals_masked are out of order.
    cals_masked = Any[nothing for _ in cals]
    Threads.@threads :dynamic for i in eachindex(cals)
        cal = cals[i]
        calmasked = deepcopy(cal)
        xs = axes(calmasked,1) .- calmasked["STAR-X"]
        ys = axes(calmasked,2) .- calmasked["STAR-Y"]
        rs = @. sqrt(xs^2 + ys'^2)
        calmasked[rs .<= psf_mask_radius] .= NaN
        cals_masked[i] = copyheader(cal, mapwindow(median, calmasked, (3,3)))
    end

    # remove bad pixels using a median filter
    skys = map(skys) do sky
        copyheader(sky,mapwindow(median, sky[:,:], (3,3)))
    end

    # We need to subract off the median values of all images for the least squares to work well.
    # Measure the median of science images from their masked counterparts (so star doesn't bias)
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
    skys = broadcast(skys, meds) do sky, med
        copyheader(sky, (sky .- med)[:,:])
    end

    if savedata && write_calmasked
        @warn "TODO: adapt file name" 
        # @progress "Writing" for (i,frame) in enumerate(cals_masked)
        #     AstroImages.writefits(@sprintf("%s.%04d.3.cal.masked.fits.gz", , i),frame)
        # end
    end


    
    @info "Calibrating raw images (sky and/or choppping)"
    savers = []
    Threads.@threads :dynamic for i in eachindex(cals)
        outfname = outfnames[i]
        A = cals[i]
        
        # If NIRC2 and L or M band, then use Jayke's BG ADI technique on the skys
        # Use two sets of BG references: sky images, and sky images transformed to match rotator position
        if occursin("NIRC2", get(header(cals[i]), "CURRINST", "")) && cals[i]["WAVE"] > 3.0e-6
            @info "Preparing extra sky BG refs to account for K-mirror position"
            use_skys = [
                skys;
                [
                    prepare_nirc2_kmirror_rotated_bgrefs(A, sky)
                    for sky in skys
                    if haskey(sky, "ROTPPOSN")
                ]
            ]
        else
            use_skys = skys
        end

        # If doing chopping, use other calibrated science images as BG references if the star is 
        # far enough away
        if chopped
            # If NIRC2 and L or M band, then use Jayke's BG ADI technique as well.
            # Use two sets of BG references: chopped science images,
            # and chopped science images transformed to match rotator position.
            if occursin("NIRC2", get(header(cals[i]), "CURRINST", "")) && cals[i]["WAVE"] > 3.0e-6
                @info "Preparing extra chopped BG refs to account for K-mirror position"
                use_cals = [cals; prepare_nirc2_kmirror_rotated_bgrefs.(Ref(A), cals)]
                use_cals_masked = [cals_masked; prepare_nirc2_kmirror_rotated_bgrefs.(Ref(A), cals_masked)]
            else
                # Otherwise proceed without the extra K mirror transformed references.
                use_cals = cals
                use_cals_masked = cals_masked
            end

            # Determine which frames are far enough away to use as a BG reference
            star_x_y = map(use_cals) do cal
                cal["STAR-X"], cal["STAR-Y"]
            end
            this_pos = star_x_y[i]
            distances = map(xy->sqrt((this_pos[1] - xy[1])^2 + (this_pos[2] - xy[2])^2), star_x_y)
            
            # These are the frames we will be allowed to use as bg references
            ii_allowed = findall(distances .>= sky_psf_radius_px)
            # Both regular sky frames and light frames that are far enough away to be used as sky backgrounds (chopping)
            chop_skys = vcat(use_cals[ii_allowed], use_skys)
            masked_chop_skys = vcat(use_cals_masked[ii_allowed], use_skys)
        else
            # If not doing chopping, just use sky images and not other calibrated science images
            chop_skys = use_skys
            masked_chop_skys = use_skys
        end
        # Rank by correlation. The user can allow us to only select N most correlated BG references
        correlation_vector = map(masked_chop_skys) do B
            mask = isfinite.(A) .& isfinite.(B)
            if count(mask) == 0 || all(A[mask] .== B[mask])
                return 0.0
            end
            return abs(cor(A[mask],B[mask],))
        end
        II = sortperm(correlation_vector, rev=true)[1:min(end,sky_N_most_similar)]
        
        # Alternative: rank most similar images purely by distance in time
        # t0 = DateTime(cal.headers["DATE-OBS"]*" "*cal.headers["UTC"],dateformat"yyyy-mm-dd HH:MM:SS.s")
        # correlation_vector = map(skysmasked) do sky
        #     ti = DateTime(sky.headers["DATE-OBS"]*" "*sky.headers["UTC"],dateformat"yyyy-mm-dd HH:MM:SS.s")
        #     return abs(t0-ti)
        # end 
        # II = sortperm(correlation_vector)[1:min(end,sky_N_most_similar)]

        # Now perform the BG calibration


        skycube = stack(chop_skys[II])
        skycube[.!isfinite.(skycube)] .= 0

        skycubemasked = stack(masked_chop_skys[II])
        skycubemasked[.!isfinite.(skycubemasked)] .= 0

        # Tile the image into large chunks, doing the fitting one chunk at a time instead 
        # of finding a global solution
        newsky = deepcopy(A)
        fill!(newsky, 0)

        tiles_x = range(1, stop=size(skycubemasked,1), length=tile_count+1)
        tiles_x_start = floor.(Int, tiles_x[begin:end-1])
        tiles_x_stop = floor.(Int, tiles_x[begin+1:end])
        tiles_y = range(1, stop=size(skycubemasked,2), length=tile_count+1)
        tiles_y_start = floor.(Int, tiles_y[begin:end-1])
        tiles_y_stop = floor.(Int, tiles_y[begin+1:end])
        j = 0
        # Go through tiles
        for (start_x,stop_x) in zip(tiles_x_start,tiles_x_stop), 
            (start_y,stop_y) in zip(tiles_y_start,tiles_y_stop)
            tilemask = (start_x:stop_x,start_y:stop_y)
            j += 1

            # speed it up by reducing the image and cube dimensions 2X in spatial coorinates
            # This is a balancing act: the more it is downsized, the more it will focus on large
            # scale structure at the expense of pixel to pixel structure.
            # skycubemasked_restrict = skycubemasked
            # cal_masked_restrict = cals_masked[i]
            skycubemasked_restrict = AstroImages.restrict(skycubemasked[tilemask...,:],(1,2))
            cal_masked_restrict = AstroImages.restrict(cals_masked[i][tilemask...],(1,2))
            
            skycubemasked_restrict[.!isfinite.(skycubemasked_restrict)] .= 0
            cal_masked_restrict[.!isfinite.(cal_masked_restrict)] .= 0

            finite_mask = isfinite.(cal_masked_restrict) .& all(isfinite, skycubemasked_restrict ,dims=3)[:,:]
            if count(finite_mask) == 0
                @warn "empty tile, no bg subtraction" tilemask
                continue
            end

            # Also ignore pixels outside a certain quantile range
            l,u = quantile(@view(cal_masked_restrict[finite_mask]), (0.001, 0.999))
            finite_mask .&= l .<= cal_masked_restrict .<= u
            l,u = quantile(@view(skycubemasked_restrict[finite_mask,:]), (0.001, 0.999))
            finite_mask .&= all(l .<= skycubemasked_restrict .<= u,dims=3)[:,:]
            finite_mask = BitMatrix(finite_mask[:,:])

            # Now calculate the coefficients of the subtraction using least-squares
            S = reshape(skycubemasked_restrict[finite_mask,:],count(finite_mask), :)
            b = vec(cal_masked_restrict[finite_mask])
            coeff = S \ b
            if any(!isfinite, coeff)
                @warn "non finite sky least square coefficients" 
                coeff .= 1 / length(coeff)
            end

            # Apply subtraction using un-masked versions
            newsky[tilemask...] = sum(skycube[tilemask...,:] .* reshape(coeff,1,1,:),dims=3)[:,:]
        end
    
        corrected = copyheader(A, A .- newsky)

        push!(corrected, History, "$(Date(Dates.now())): LOCI background subtraction.")
        corrected["BGTILEN"] = tile_count


        # push!(savers, (outfname, corrected))
        push!(savers, (outfname, corrected))
        # AstroImages.writefits(outfname,corrected)
        println(outfname, "\t($i) calc")
    end

    for (fname, img) in savers
        AstroImages.writefits(fname,img)
        println(fname, "\t write")
    end


    return


end
export bgsub2