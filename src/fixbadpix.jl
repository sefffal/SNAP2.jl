export fixbadpix

# TODO: we need this to actually work off the derivative change, not abs deviation.
function fixbadpix(fnames_pattern::AbstractString,force=false)

    fnames = glob(fnames_pattern)
    imgs = load.(fnames)
    i = 0
    for (fname, targ) in zip(fnames, imgs)
        i+=1
        outfname = replace(fname, ".fits"=>".bpfix.fits")
        if !force && isfile(outfname) &&  Base.Filesystem.mtime(outfname) > Base.Filesystem.mtime(fname)
            continue
        end
        out = fixbadpix(targ)
        AstroImages.writefits(outfname, out)
    end

end
function fixbadpix(img::AstroImage)
    l = 3
    for I in CartesianIndices(img)
        ix = I[1]
        iy = I[2]
        # Get patch
        patch = img[max(begin,ix-l):min(end,ix+l),max(begin,iy-l):min(end,iy+l),]
        # remove this pixel
        px  = img[ix,iy]
        other_px = setdiff(vec(patch), px)
        if abs(px) > 5std(filter(isfinite,other_px))
            img[ix,iy] = NaN
        end
    end
    return img
end
