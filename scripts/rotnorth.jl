#=
This script prepares frames for subtraction.

One script rotates images North-up.

Another script rotates and/or scales reference frames to match a given target frame.

=#

using ImageTransformations
using CoordinateTransformations: CoordinateTransformations
using Rotations: Rotations
using LinearAlgebra
function rotnorth(
    pattern::AbstractString,
    offset=[0,0];
    force=false,
)
    fnames = Glob.glob(pattern)
    outfnames = String[]
    for fname in fnames
        outfname = replace(fname,".fits"=>".rotnorth.fits")
        if !force && isfile(outfname) && Base.Filesystem.mtime(outfname) > Base.Filesystem.mtime(fname)
            push!(outfnames, outfname)
            continue
        end
        img = load(fname)
        cx, cy = (img["STAR-X"], img["STAR-Y"]) .+ offset

        # TODO: must confirm the rotation direction is correct!
        rot_deg = -img["ANGLE_MEAN"]
        rot_rad = deg2rad(rot_deg)

        tfrm = CoordinateTransformations.LinearMap(
            Rotations.RotMatrix(rot_rad)
        )
        tfrm = CoordinateTransformations.recenter(tfrm, (cx,cy))
        applied_dat = collect(ImageTransformations.warp(
            collect(img),
            tfrm,
            axes(img),
            # Broken for some reason:...
            # Interpolations.Quadratic(Interpolations.Flat(Interpolations.OnGrid())),
            NaN,
        ))

        applied = copyheader(img, applied_dat)
        push!(applied, History, "$(Date(Dates.now())): Rotated North-Up.")
        applied["ORIG.ANGLE_MEAN"] = applied["ANGLE_MEAN"]
        applied["ANGLE_MEAN"] = 0
        AstroImages.writefits(outfname, applied)
        println(outfname, "\t rotnorth ", rot_deg)
        push!(outfnames, outfname)
    end
    return outfnames
end
