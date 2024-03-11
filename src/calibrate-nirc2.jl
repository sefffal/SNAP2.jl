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

export calibrate_nirc2

"""
    calibrate("sequences/abc/profile.toml")

Perform NIRC2 calibration. 
see `calibrate(conf_file)` for public interface.

Returns named tuple with:
 - `registered`, an array of calibrated & registered images.
"""
function calibrate_nirc2(conf_fname; verbose=true, savedata=true, showplots=true, force=false)

    verbose && println(SNAP.banner()*"\n\nCALIBRATE")
   
    # Read the configuration
    conf = TOML.parsefile(conf_fname)
    cal = conf["calibrate"]

    # Set up some nice file paths 
    # Create a nice sanitized name for reporting these results
    if !occursin('/', conf_fname)
        conf_fname = joinpath(".", conf_fname)
    end
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

    ## Load calibration files.
    # Can be a list of paths, or a shell glob pattern.
    function loadfiles(paths::AbstractArray)
        if isempty(paths)
            @warn "No files found with paths $paths"
        end
        @progress "  loading" files = [
            SNAP.header_normalize!(load(path))
            for path in paths
        ]
        return files
    end
    function loadheaders(paths::AbstractArray)
        map(paths) do path
            FITS(path) do fit
                read_header(first(fit))
            end
        end
    end


  
    loadfiles(pathspec::AbstractString) = loadfiles(glob(pathspec))

    # Load files and linearize them using the above routine
    function loadfiles_lin(pathspec::Union{<:AbstractString,AbstractArray})
        files = loadfiles(pathspec)
        
        return SNAP.linearize_nirc2_data!!.(files)
    end

    loadheaders(pathspec::AbstractString) = loadfiles(glob(pathspec))
    function loadheader(path::AbstractString)
        FITS(path) do fit
            read_header(first(fit))
        end
    end



    ##
    function stackframes(frames, method)
        if length(frames) == 0
            error("No matching frames")
        elseif length(frames) == 1
            return only(frames)
        end
        combined = method(stack(frames),dims=3)
        return copyheader(first(frames), combined)
    end

    # TODO: Should have a sanity check to detect duplicates

    write_masters = get(cal, "write_masters", true)

    verbose && @info "Loading and stacking masters"

    ## Search for matching files
    fnames = conf["calibrate"]["raw_path"]
    fnames_dark = loadfiles_lin(cal["dark"])
    fnames_dark_psf = loadfiles_lin(cal["dark_psf"])
    fnames_flaton = loadfiles_lin(cal["flaton"])
    fnames_flatoff = loadfiles_lin(cal["flatoff"])
    fnames_sky = cal["sky"]
    fnames_sky_psf = loadfiles_lin(cal["sky_psf"])

    length(fnames_dark) == 0 && error("No matching dark frames")
    dark    = stackframes(fnames_dark, median)

    fs1 = loadfiles(fnames[1:1])[1]
    if size(dark) != size(fs1)
        @info "Cropping dark to match first data file dimensions."
        dark = crop_A_to_B(dark, fs1)
    end

    # Create bad pixel map
    bpm =float.(deepcopy(dark))
    # Can also load in a manual bad pixel map (but auto is always also applied)
    if haskey(cal, "badpix")
        verbose && @info "Loading additional manual bad pixel maps"
        # bpm.= 0
        @show cal["badpix"]
        bpm = stackframes(loadfiles(cal["badpix"]), sum)
        med = median(filter(isfinite, bpm))
        bpm = float.(crop_A_to_B(bpm, fs1))
        bpm[bpm .> 2med] .= NaN
    else
        SNAP.excludebadpixels!(bpm, 2.5, 2)
    end
    if size(bpm) != size(fs1)
        @info "Cropping bpm to match first data file dimensions."
        bpm = crop_A_to_B(bpm, fs1)
    end


    ##
    darkpsf = stackframes(fnames_dark_psf, median)[:,:]



    psfs     = loadfiles_lin(cal["raw_psf_path"])

    # Since dark exposure does not necessarily match the other exposures,
    # ideally we want to do dark scaling. But mostly we do not have a bias.
    # Instead, the dark is effectively the bias frame and the sky is the "dark"


        
    length(fnames_sky_psf) == 0 && @warn("No matching sky_psf frames")
    if length(fnames_sky_psf) > 0
        skypsfraw = stackframes(fnames_sky_psf, median)
        # We use the PSF skies as flats
        skypsf = copyheader(
            skypsfraw,
            skypsfraw .- darkpsf,
        )
    else
        skypsfraw = deepcopy(darkpsf)
        fill!(skypsfraw, 1)
        # We use the PSF skies as flats
        skypsf = skypsfraw
    end
    flatmedpsf = median(filter(isfinite, skypsfraw))
    


    ###
    verbose && @info "Calibrating PSF images"
    bpm_psf = SNAP.crop_A_to_B(bpm, first(psfs))
    # calibrated_psfs = parmap(psfs) do psf
    calibrated_psfs = map(psfs) do psf
        dat = (psf .- darkpsf) ./
            # psf.headers["TOTEXP"]
            skypsf .* flatmedpsf .+
            bpm_psf
        dat = dat[:,:,1,]
        # Remove any weird bias offset between PSF and data
        dat = dat .- median(filter(isfinite,[collect(dat[begin:begin+3,:]) collect(dat[end-3:end,:])]))
        psfcal = copyheader(psf,dat)
        psfcal = SNAP.interpbadpixels!(psfcal)
        psfcal[.!isfinite.(psfcal)] .= median(filter(isfinite, psfcal))
        return psfcal
    end
    if write_masters && savedata
        for i in eachindex(calibrated_psfs) 
            try
                AstroImages.writefits("$caldir.psf.cal.$i.fits", calibrated_psfs[i])
            catch err
                @error "Error writing file:" exception=err
                rethrow(err)
            end
        end
    end
    verbose && @info "Transforming PSFs"
    output_axes_psf = (1:80, 1:80)
    star_x = round(Int, get(conf["calibrate"], "star_x", 0))
    star_y = round(Int, get(conf["calibrate"], "star_y", 0))
    transformed_psfs = map(calibrated_psfs) do psf
        if star_x != 0
            results = SNAP.centrefit_fixedstd(psf, (star_x-50, star_x+50), (star_y-50, star_y+50), std=3)
        else
            results = SNAP.centrefit_fixedstd(psf, (1, size(psf,1)), (1, size(psf,2)), std=3)
        end
        tfrm = CoordinateTransformations.Translation((results.x,results.y).- mean.(output_axes_psf))
        applied = parent(ImageTransformations.warp(
            collect(psf),
            tfrm,
            output_axes_psf,
            Interpolations.Quadratic(Interpolations.Flat(Interpolations.OnGrid())),
            NaN,
        ))
        return copyheader(psf, applied)
    end
    if write_masters && savedata
        for i in eachindex(transformed_psfs) 
            AstroImages.writefits("$caldir.psf.reg.$i.fits", transformed_psfs[i])
        end
    end
    verbose &&  @info "Stacking PSF"
    # psf = stackframes(transformed_psfs, mean)
    psf_data = zeros(eltype(first(transformed_psfs)),size(first(transformed_psfs)))
    for I in eachindex(psf_data)
        pixels = filter(isfinite, [psf[I] for psf in transformed_psfs])
        if length(pixels) > 0
            psf_data[I] = mean(pixels)
        end
    end
    psf = copyheader(first(transformed_psfs), psf_data)
    showplots && display(imview(psf))

    # Create a copy without the central core for using as a cross correlation with saturated images
    psf_cored = deepcopy(psf)
    r = imgsep(psf_cored)
    psf_cored[r .< 8] .= 0

    savedata && verbose && @info "Saving PSF"
    savedata && AstroImages.writefits(cal["output_psf_path"],psf)
    savedata && verbose && @info "Saving cored PSF"
    savedata && AstroImages.writefits(replace(cal["output_psf_path"], ".fits"=>"-cored.fits"), psf_cored) 


    # If a pixel exceeds 4sigma the surrounding 6x6 pixel box (2 pixels on either side),
    # then mark it as bad.
    bpm[isfinite.(bpm)] .= 0.0

    dark = SNAP.interpbadpixels!(dark .+ bpm)

    verbose && @info "Compiling flat frames"
    fnames_flaton = map(fnames_flaton) do img
        img = crop_A_to_B(img, fs1)
        # Bad pixel map has 0 everywhere except NaN on bad pixels
        img = img.+bpm
        SNAP.interpbadpixels!(img)
        return img
    end
    fnames_flatoff = map(fnames_flatoff) do img
        img = crop_A_to_B(img, fs1)
        # Bad pixel map has 0 everywhere except NaN on bad pixels
        img = img.+bpm
        SNAP.interpbadpixels!(img)
        return img
    end

    length(fnames_dark_psf) == 0 && error("No matching dark_psf frames")
    length(fnames_flaton) == 0 && error("No matching flaton frames")
    length(fnames_flatoff) == 0 && error("No matching flatoff frames")



    flaton  = stackframes(fnames_flaton, median)[:,:]
    flatoff = stackframes(fnames_flatoff, median)[:,:]


    # Flat is difference between flat lamp on and off.
    # No need to dark subtract
    flat = copyheader(flaton, flaton .- flatoff)
    
    # Calibrate sky using dark and flat
    flatmed = median(filter(isfinite, flat))

    if any(==(0), flat)
        @warn "zero values present in flats! This may lead to poor results"
        flat[flat .== 0] .= flatmed
    end

        
    if get(cal, "skipsky", false) || length(fnames_sky) == 0
        sky = deepcopy(dark)
        fill!(sky, 0)
        sky = [sky]
        @warn("No matching sky frames or skipsky=true")
    else
        # We support grouped batches of sky filenames stacked together, and then LOCI subtracted
        if eltype(fnames_sky) <: AbstractString
            skys = [[sky] for sky in loadfiles_lin(fnames_sky)]
        else
            skys = map(fnames_sky) do fnames_group
                loadfiles_lin(fnames_group)
            end
        end
        skyraw = map(fnames->stackframes(fnames, median), skys)
        @info "Sky groups" N=length(skyraw)
        sky = map(
                skyraw->copyheader(
                    skyraw,
                    (skyraw .- dark) ./ flat .* flatmed .+ bpm,
                ),
                skyraw
        )
        sky = map(SNAP.interpbadpixels!,sky)
        # sky.data ./= sky.headers["TOTEXP"]
    end
    if size(flat) != size(fs1)
        @info "Cropping flat to match first data file dimensions."
        flat = crop_A_to_B(flat, fs1)
    end
    sky = map(sky) do s
        if size(s) != size(fs1)
            @info "Cropping sky frame to match first data file dimensions."
            s = crop_A_to_B(s, fs1)
        end
        return s
    end


    ## Optionally write out copies of the masters
    if write_masters && savedata
        verbose && @info "Writing masters"
        AstroImages.writefits("$caldir.dark.fits.gz", dark)
        AstroImages.writefits("$caldir.darkpsf.fits.gz", darkpsf)
        AstroImages.writefits("$caldir.skypsf.fits.gz", skypsf)
        AstroImages.writefits("$caldir.flat.fits.gz", flat)
        AstroImages.writefits("$caldir.sky.fits.gz", sky...)
        AstroImages.writefits("$caldir.bpm.fits.gz", bpm)
    end


    verbose && @info "Reading files"
    # Recommend you test with one or two frames to begin
    fnames = cal["raw_path"]


    fnames_out = map(fnames) do fname
        # Move files to cal subdirectory, and append "cal" to filename
        fname = caldir*"."*splitpath(fname)[end]
        fname = replace(fname, ".fits"=>".cal.fits")
        # We'll want to use compressed files going forward.
        # We often pad out the images with NaN eg after registration
        # and these take up tons of space (but almost none when compressd)
        if !endswith(fname, ".gz")
            fname = fname*".gz"
        end
        return fname
    end
    fnames_II = filter(eachindex(fnames)) do i
        if !force && isfile(fnames_out[i])  && Base.Filesystem.mtime(fnames_out[i]) > Base.Filesystem.mtime(fnames[i])
            return false
        end
        return true
    end
    if isempty(fnames_II)
        return
    end

    # TODO: better to go through the frames one at a time, loading, calibrating, then
    # writing out.
    # Need to queue up loading and saving on a single thread though.

    ### WARNING: FITSIO not THREADSAFE?? ##
    # Load the files first, then process them in parallel.
    raws = loadfiles_lin(fnames[fnames_II])

    for raw in raws
        SNAP.interpbadpixels!(raw .+ bpm)
    end

    # if savedata && get(conf["calibrate"], "write_calibrated", true)
    #     verbose && @info "Writing bad-pixel corrected files for second bg subtraction"
    #     @progress "Writing" for (i,frame) in zip(fnames_II, raws)
    #         AstroImages.writefits(@sprintf("%s.%04d.1.badpixfix.fits.gz", caldir, i),frame[:,:])
    #     end
    # end
    # gif_fname = "$caldir.anim.1.badpixfix.gif"
    # @ANIMATE
    # frames = stack(raws)[begin:4:end,begin:4:end,:]
    # save(gif_fname, imview(frames), fps=10)

    verbose && @info "Calibrating raw images (dark, flat, bad-pix, not sky)"
    to_save = []
    # Threads.@threads :dynamic 
    for (i,raw) in collect(zip(fnames_II, raws))
        data = (@. (raw - dark) / flat * flatmed .+ bpm)[:,:]
        # cal = RasterImage(raw; data=data, xcoords, ycoords)
        cal = copyheader(raw, data)
        push!(cal, History, "$(Date(Dates.now())): dark & flat calibrated")

        # This line is the bottle-neck by factor 100
        cal = SNAP.interpbadpixels!(cal)
        push!(cal, History, "$(Date(Dates.now())): replaced bad pixels with neighbours")

        # Apply the distortion solution
        dat = SNAP.nirc2_dewarp!(AstroImage(collect(cal)))
        cal = copyheader(cal, dat)

        verbose && println(fnames_out[i], "\t", cal["DATE-OBS"], " ", cal["UTC"], "\t($i)")
        push!(to_save, (fnames_out[i],cal[:,:]))
        # saver = @async AstroImages.writefits(fnames_out[i],cal[:,:])
        # push!(savers,saver)
    end
    @info "Writing"
    for (f,d) in to_save
        AstroImages.writefits(f,d)
    end
    # fetch.(savers)

    # gif_fname = "$caldir.anim.2.initial-cal.gif"
    # @ANIMATE
    # frames = stack(calibrated[begin:4:end,begin:4:end])
    # save(gif_fname, imview(frames), fps=10)

    return
end




"""
ridge_regression(Y,X,λ,β₀=0)

Calculate ridge regression estimate with target vector β₀.

https://discourse.julialang.org/t/hello-i-am-stuck-with-a-problem-regarding-linear-regression-and-ridge-regression/75382/11?u=sefffal
"""
function ridge_regression(Y,X,λ,β₀=0)
    K = size(X,2)
    isa(β₀,Number) && (β₀=fill(β₀,K))
    local b
    try
        b = (X'X+λ*I)\(X'Y+λ*β₀)      #same as inv(X'X+λ*I)*(X'Y+λ*β₀)
    catch err
        @error "error in ridge ridge_regression" exception=(err, catch_backtrace())
        return NaN .* (X \ Y)
    end
    return b
end


