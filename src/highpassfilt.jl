export highpassfilt
using ImageFiltering
function highpassfilt(img::AstroImageMat; k=Kernel.gaussian(10))
    mask_finite = isfinite.(img)
    temp = deepcopy(img)
    temp[.!mask_finite] .= 0
    out = temp .- imfilter(mapwindow(median,temp,(3,3)), k)
    out[.!mask_finite] .= NaN
    return out
end

function highpassfilt(pattern::AbstractString; force=false, k=Kernel.gaussian(10))
    fnames = Glob.glob(pattern)
    return highpassfilt(fnames; force, k)
end
function highpassfilt(fnames::AbstractVector{<:AbstractString}; force=false, k=Kernel.gaussian(10))
    for fname in fnames
        # Calculate contrast if not there, or load
        outfname = replace(fname, ".fits"=>".hpf.fits")
        # TODO: check itime
        if !force && isfile(outfname) && Base.Filesystem.mtime(outfname) > Base.Filesystem.mtime(fname)
            continue
        end
        img = load(fname)
        out = highpassfilt(img; k)
        AstroImages.writefits(outfname, copyheader(out, collect(out)))
        println(outfname)
    end
end