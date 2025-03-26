using Optim
using AstroImages
using ImageFiltering
using CoordinateTransformations, ImageTransformations
using ImageTransformations.Interpolations
using FITSIO

export register

function register(pattern::AbstractString, template_fname::AbstractString; kwargs...)
    fnames = Glob.glob(pattern)
    return register(fnames, template_fname; kwargs...)
end
function register(fnames::AbstractArray{<:AbstractString}, template_fname::AbstractString; force=false, cropsize=90, searchsize=cropsize)
    template = load(template_fname)

    # high pass filter and median filter the template
    newtemplate = mapwindow(median,template,(3,3)) .- imfilter(template, Kernel.gaussian(5))
    # newtemplate .-= minimum(newtemplate)
    # newtemplate = sqrt.(newtemplate)
    # newtemplate = imfilter(newtemplate, Kernel.gaussian(1.0))
    

    # Output axes should be large enough to cover all edges of the images given their star positions.
    # Make it an odd number so centre pixel is star location.
    outsize = maximum(map(fnames) do fname
        FITS(fname) do fits
            sz = size(fits[1])[1:2]
            h = read_header(fits[1])
            xy = h["STAR-X"], h["STAR-Y"]
            # need distances to both corner.
            topleft_dist = norm(xy)
            topright_dist = norm(sz .- xy)
            ceil(Int, 2√2 * max(topleft_dist, topright_dist))
        end
    end)
    
    # We want the output size to be identical for "similar" reductions,
    # so find the nearest power of two that meets this criteria.
    # This is just a heuristic that makes it more convenient to compare 
    # results from different sequences, or results from slightly different
    # registration results for the same sequence.
    powtwo = round.(Int, 2 .^ (5:0.5:11))
    i = findfirst(>(0), powtwo .- outsize)
    if isnothing(i)
        i = lastindex(powtwo)
    end
    outsize = powtwo[i]+1
    ax_out = (1:outsize,1:outsize)
    @info "Determined output canvas size" ax_out


    frames = AstroImage[]
    fnames_out = String[]
    @info "Loading frames"
    for fname in fnames
        outfname = replace(fname, ".fits"=>".reg.fits",)
        if force || !isfile(outfname) || Base.Filesystem.mtime(fname) > Base.Filesystem.mtime(outfname)  || Base.Filesystem.mtime(template_fname) > Base.Filesystem.mtime(outfname)
            push!(frames, load(fname))
            push!(fnames_out, outfname)
        end
    end

    savers = []
    Threads.@threads :dynamic  for (satimg, outfname) in collect(zip(frames,fnames_out))
        x_guess = round(Int, satimg["STAR-X"])
        y_guess = round(Int, satimg["STAR-Y"])
        crop_left = max(first(axes(satimg,1)),x_guess-cropsize)
        crop_right = min(last(axes(satimg,2)), x_guess+cropsize)
        crop_bottom = max(first(axes(satimg,1)),y_guess-cropsize)
        crop_top = min(last(axes(satimg,2)),y_guess+cropsize) # 👕
        c = satimg[crop_left:crop_right, crop_bottom:crop_top]

        # High pass filter
        cmask = mapwindow(median,c,(3,3)) .- imfilter(c, Kernel.gaussian(5))
        # cmask .-= minimum(cmask)
        # cmask = sqrt.(cmask)
        # cmask = imfilter(cmask, Kernel.gaussian(1.0))
    
        q = imfilter((cmask), ImageFiltering.centered((newtemplate)))
        #display(imview(q))

        std = 4
        gauss2D(A, offset, μx, μy, σx, σy, x, y) =
            A * exp(-(((x - μx)^2) / (2 * σx^2) + ((y - μy)^2) / (2 * σy^2))) + offset

        model(A, offset, μx, μy, x, y) = gauss2D(A, offset, μx, μy, std, std, x, y)

        function objective(params)
            A, offset, μx, μy = params
            result = 0.0
            # try # sig might be infeasible so we have to handle this case
                for I in CartesianIndices(q)
                    x, y = Tuple(I)
                    est = model(A, offset, μx, μy, x, y)
                    meas = q[I]
                    resid = (meas - est)^2
                    # Skip NaN pixels in image - they give no weight.
                    if isfinite(resid)
                        result += resid
                    end
                end
            # catch e
                # If the gaussian throws an error.
                # result = Inf
            # end
            return log(result)
        end

        startI = argmax(q)
        guess_offset = mean(q)
        guess = [
            # A, offset, μx, μy
            # maximum(pixels) - mean(pixels),
            q[startI],
            guess_offset,# mean(pixels),
            x_guess-crop_left,
            y_guess-crop_bottom,
        ]
        lower = [0, -1e32, x_guess-crop_left-searchsize, y_guess-crop_bottom-searchsize, ]
        upper = [1e32, 1e32, x_guess-crop_left+searchsize, y_guess-crop_bottom+searchsize, ]
        # test for an example starting point
        result = optimize(
            objective,
            lower,
            upper,
            guess,
            Optim.Fminbox(Optim.LBFGS()),
            Optim.Options(
                show_trace=false,
                time_limit=10.0, # Soft upper limit in seconds before returning best guess.
                allow_f_increases=true
            ),
            autodiff=:forward,
        )
        if !Optim.converged(result)
            @error("Centre-fit did not converge")
            display(result)
        end

        # Now that we have found the centre of our cropped image, we can apply the full translation to
        # 

        
        tfrm = CoordinateTransformations.Translation((result.minimizer[3]+x_guess-cropsize-0.5,result.minimizer[4]+y_guess-cropsize-0.5))
        applied = collect(ImageTransformations.warp(
            collect(satimg)[:,:],
            tfrm,
            (-outsize÷2:outsize÷2, -outsize÷2:outsize÷2),
            # Broken for some reason:...
            # Interpolations.Quadratic(Interpolations.Flat(Interpolations.OnGrid())),
            NaN,
        ))

        registered = copyheader(satimg, applied)
        push!(registered, History, "$(Date(Dates.now())): Registered star to image centre.")
        registered["STAR-X"] = mean(ax_out[1])
        registered["STAR-Y"] = mean(ax_out[2])
        registered["CRPIX1"] = mean(ax_out[1])
        registered["CRPIX2"] = mean(ax_out[2])
        
        push!(savers, (outfname, registered))
        println(outfname, "\t register ", tfrm)
    end
    @info "Writing files"
    for (fname, img) in savers
        AstroImages.writefits(fname, img)
    end

    return
end


"""
Version of `register` that just uses the saved "STAR-X" and "STAR-Y"
header keywords to re-position the images (doesn't do a cross correlation
with a template and a gaussian fit).
"""
function register_nofit(pattern; force=false)
    fnames = Glob.glob(pattern)

    # Output axes should be large enough to cover all edges of the images given their star positions.
    # Make it an odd number to centre pixel is star location.
    outsize = maximum(map(fnames) do fname
        FITS(fname) do fits
            sz = size(fits[1])[1:2]
            h = read_header(fits[1])
            xy = h["STAR-X"], h["STAR-Y"]
            # need distances to both corner.
            topleft_dist = norm(xy)
            topright_dist = norm(sz .- xy)
            ceil(Int, 2√2 * max(topleft_dist, topright_dist))
        end
    end)
    
    # We want the output size to be identical for "similar" reductions,
    # so find the nearest power of two that meets this criteria.
    powtwo = 2 .^ (5:11)
    i = findfirst(>(0), powtwo .- outsize)
    if isnothing(i)
        i = lastindex(powtwo)
    end
    outsize = powtwo[i]+1
    ax_out = (1:outsize,1:outsize)

    frames = AstroImage[]
    fnames_out = String[]
    @info "Loading frames"
    for fname in fnames
        outfname = replace(fname, ".fits"=>".reg.fits",)
        if force || !isfile(outfname) || Base.Filesystem.mtime(fname) > Base.Filesystem.mtime(outfname)
            push!(frames, load(fname))
            push!(fnames_out, outfname)
        end
    end

    savers = []
    Threads.@threads :dynamic for (satimg, outfname) in collect(zip(frames,fnames_out))
        x_guess = round(Int, satimg["STAR-X"])
        y_guess = round(Int, satimg["STAR-Y"])
        
        tfrm = CoordinateTransformations.Translation((x_guess-0.5,y_guess-0.5))
        applied = collect(ImageTransformations.warp(
            collect(satimg)[:,:],
            tfrm,
            (-outsize÷2:outsize÷2, -outsize÷2:outsize÷2),
            # Broken for some reason:...
            # Interpolations.Quadratic(Interpolations.Flat(Interpolations.OnGrid())),
            NaN,
        ))

        registered = copyheader(satimg, applied)
        push!(registered, History, "$(Date(Dates.now())): Registered star to image centre (used STAR-X/Y not XCORR).")
        registered["STAR-X"] = mean(ax_out[1])
        registered["STAR-Y"] = mean(ax_out[2])
        registered["CRPIX1"] = mean(ax_out[1])
        registered["CRPIX2"] = mean(ax_out[2])
        
        push!(savers, (outfname, registered))
        println(outfname, "\t register ", tfrm)
    end
    @info "Writing files"
    for (fname, img) in savers
        AstroImages.writefits(fname, img)
    end

    return
end