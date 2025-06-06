using Statistics, StatsBase
using Interpolations

export stackframes, stackframes_contweight

"""
    stackframes_contweight(func, filepaths)
    stackframes_contweight(filepaths)

Given a glob pattern matching one or more files, stack the files weighted by contrast.
Function defaults to a contrast-weighted median.
"""
function stackframes_contweight(pattern::AbstractString; force=false,step=4)
    fnames = Glob.glob(pattern)
    outname = replace(pattern, "*"=>"_", ".fits"=>".stackedw.fits")
    return stackframes_contweight(median, fnames, outname; force, step)
end
function stackframes_contweight(method::Base.Callable, pattern::AbstractString; force=false,step=4)
    fnames = Glob.glob(pattern)
    outname = replace(pattern, "*"=>"_", ".fits"=>".stackedw.fits")
    return stackframes_contweight(method, fnames, outname; force, step)
end
function stackframes_contweight(method::Base.Callable, fnames::Vector{<:AbstractString}, outname::AbstractString;force=false,step=4)

    # Check if there is any work to do
    if !force && isfile(outname) && Base.Filesystem.mtime(outname) > maximum(Base.Filesystem.mtime.(fnames))
        out = load(outname)
        return out
    end
    imgs = convert(Vector{AstroImage}, load.(fnames))

    contrasts = contrast(fnames; force, step)

    return stackframes_contweight(method, imgs, contrasts, outname)
end
function stackframes_contweight(method, imgs::Vector{<:AstroImage}, contrasts::Vector, outname::AbstractString)
    # Could think about ways to do this with OnlineStats so we can do rolling
    # medians, means in constant time with new frames

    contrast_interpolators = map(contrasts) do (sep,cont)
        return Interpolations.linear_interpolation(
            sep, log.(cont),
            extrapolation_bc=Line()#convert(eltype(cont), NaN)
        )
    end

    img_seps = imgsep.(imgs)

    # Place into 3D cube
    cube = stack(imgs)
    out = similar(first(imgs))

    # Use a function barrier here for performance, otherwise `imgs` is Any.
    (function (out, img_seps, cube, contrast_interpolators)
        fill!(out, NaN)
        mask_dat = falses(size(cube,3))
        mask_cont = falses(size(cube,3))
        conts = zeros(eltype(cube),size(cube,3))
        for I1 in axes(out,1), I2 in axes(out,2)
            dat = @view cube[I1,I2,:]
            mask_dat .= isfinite.(dat)
            if count(mask_dat) == 0
                continue
            end
            for I3 in axes(cube,3)
                conts[I3] = 1 / exp(contrast_interpolators[I3](img_seps[I3][I1,I2]))^2
            end
            mask_cont .= isfinite.(conts)
            # If all data has non-finite contrast, weight equally instead of ignoring.
            if count(mask_cont) == 0
                continue
                # conts .= 1
                # mask_cont .= true
            end
            mask_dat .&= mask_cont
            if count(mask_dat) == 0
                continue
            end
            weights = AnalyticWeights(view(conts,mask_dat))
            out[I1,I2] = method(view(dat,mask_dat), weights)
        end
        return out
    end)(out, img_seps, cube, contrast_interpolators)

    push!(out, History, "$(Date(Dates.now())): Stacked frames (contrast weighted).")
    AstroImages.writefits(outname, out)
    println(outname)
    return out
end

"""
    stackframes(func, filepaths)
    stackframes(filepaths)

Given a glob pattern matching one or more files, stack the files weighted by contrast.
Function defaults to a contrast-weighted median.
"""
function stackframes(pattern::AbstractString; force=false)
    return stackframes(median, pattern; force)
end
function stackframes(method::Base.Callable, pattern::AbstractString; force=false)
    fnames = Glob.glob(pattern)
    outname = replace(pattern, "*"=>"_", ".fits"=>".stacked.fits")
    return stackframes(method, fnames, outname; force)
end
function stackframes(method::Base.Callable, fnames::AbstractVector{<:AbstractString}, outname::AbstractString; force=false)

    # Could think about ways to do this with OnlineStats so we can do rolling
    # medians, means in constant time with new frames
    # Check if there is any work to do
    if !force && isfile(outname) && Base.Filesystem.mtime(outname) > maximum(Base.Filesystem.mtime.(fnames))
        out = load(outname)
        return out
    end
    imgs = convert(Vector{AstroImageMat}, load.(fnames))
    return stackframes(method, imgs, outname; force)
end

function stackframes(method::Base.Callable, imgs::Vector{<:AstroImage}, outname::AbstractString; force=false)

    # Place into 3D cube
    cube = stack(imgs)
    out = similar(first(imgs))

    # Use a function barrier here for performance, otherwise `imgs` is Any.
    (function (out, cube)
        fill!(out, NaN)
        mask_dat = falses(size(cube,3))
        mask_cont = falses(size(cube,3))
        conts = zeros(eltype(cube),size(cube,3))
        for I1 in axes(out,1), I2 in axes(out,2)
            dat = @view cube[I1,I2,:]
            mask_dat .= isfinite.(dat)
            if count(mask_dat) == 0
                continue
            end
            out[I1,I2] = method(view(dat,mask_dat))
        end
        return out
    end)(out, cube)


    push!(out, History, "$(Date(Dates.now())): Stacked frames (un-weighted).")
    if method == median
        out["STACKFUN"] = "median"
    elseif method == mean
        out["STACKFUN"] = "mean"
    else
        out["STACKFUN"] = string(method)# will look kind of wonky for an arbitrary function, but better than nothing
    end
    out["STACKFUN", Comment] = "function used to stack data, in order to produce this image"
    AstroImages.writefits(outname, out)
    println(outname)
    return out
end