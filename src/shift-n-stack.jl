
export shiftnstack

function shiftnstack(pattern::AbstractString,xs=(-2:2),ys=(-2:2))
    fnames = glob(pattern)
    fnames_rot = replace.(fnames, (".fits"=>".rotnorth.fits",))
    outs = broadcast(xs, ys') do x, y
        println(x, " ", y)
        rotnorth(pattern, [x,y], force=true)
        fnameout = replace(pattern, "*"=>"_", ".fits"=>".shiftnstack_$(x)_$(y).fits")
        out = stackframes(median, fnames_rot, fnameout, force=true)
        return out
    end
    return outs
end
# ##
# AstroImages.writefits("cal/aflep-1-v2.*.cal.bgsub.reg.fluxnorm.sub.rotnorth.shift-n-stack.fits.gz",stack([a.out for a in outs2]))



# ##
# fnames = Glob.glob("cal/aflep-1-v2.*.cal.fits.gz")
# fs = load.(fnames)
# ang1 = getindex.(fs, "ANGLE_MEAN")

# ##
# fnames = Glob.glob("cal/aflep-1-v2.*.cal.bgsub.fits.gz")
# fs = load.(fnames)
# ang1 = getindex.(fs, "ANGLE_MEAN")
# fig = Makie.lines(ang1)

# ##
# fnames = Glob.glob("cal/aflep-1-v2.*.cal.bgsub.reg.fluxnorm.sub.fits.gz")
# fs = load.(fnames)
# ang2 = getindex.(fs, "ANGLE_MEAN")
# fig = Makie.lines(ang2)

# ##
# fs2 = SNAP.header_normalize_nirc2!.(fs)
# ang3 = getindex.(fs2, "ANGLE_MEAN")
# Makie.lines(ang3)

# ##
# fig = Makie.lines(ang1)
# Makie.lines!(ang2)
# fig