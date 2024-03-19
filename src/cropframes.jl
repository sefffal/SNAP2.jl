export cropframes
# Crop frames to +- croplen pixels from the star.
function cropframes(
    pattern::AbstractString,
    croplen::Int;
    force=false
)
    fnames = Glob.glob(pattern)
    for fname in fnames
        outfname = replace(fname,".fits"=>".crop.fits")
        if !force && isfile(outfname) && Base.Filesystem.mtime(outfname) > Base.Filesystem.mtime(fname)
            continue
        end
        img = load(fname)

        if haskey(img, "STAR-X")
            cx = round(Int,img["STAR-X"])
            cy = round(Int,img["STAR-Y"])
        else
            cx = round(Int,mean(axes(img,1)))
            cy = round(Int,mean(axes(img,2)))
        end
        out = img[cx-croplen:cx+croplen,cy-croplen:cy+croplen]
        out["STAR-X"] = croplen + 1 + cx - round(cx)
        out["STAR-Y"] = croplen + 1 + cy - round(cy)

        push!(out, History, "$(Date(Dates.now())): Cropped image to central region.")
        AstroImages.writefits(outfname, out)
        println(outfname, "\t crop ± $croplen ")
    end
end