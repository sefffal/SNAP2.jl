export binframes
using ImageFiltering
function binframes(img::AstroImage;)
    dat = restrict(collect(img), (1,2))
    out = copyheader(img, dat)
    if haskey(out, "STAR-X")
        out["STAR-X"] *= size(out,1)/size(img,1)
        out["STAR-Y"] *= size(out,2)/size(img,2)
    end
    return out
end
function binframes(pattern::AbstractString; force=false,)
    fnames = Glob.glob(pattern)
    return binframes(fnames; force,)
end
function binframes(fnames::AbstractVector{<:AbstractString}; force=false,)
    for fname in fnames
        # Calculate contrast if not there, or load
        outfname = replace(fname, ".fits"=>".bin.fits")
        # TODO: check itime
        if !force && isfile(outfname) && Base.Filesystem.mtime(outfname) > Base.Filesystem.mtime(fname)
            continue
        end
        img = load(fname)
        out = binframes(img;)
        AstroImages.writefits(outfname, copyheader(out, collect(out)))
        println(outfname)
    end
end