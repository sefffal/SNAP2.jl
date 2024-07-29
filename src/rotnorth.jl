#=
This script prepares frames for subtraction.

One script rotates images North-up.

Another script rotates and/or scales reference frames to match a given target frame.

=#

using ImageTransformations
using CoordinateTransformations: CoordinateTransformations
using Rotations: Rotations
using LinearAlgebra

export rotnorth, rotrestore

function rotnorth(
    pattern::AbstractString,
    offset=[0,0];
    kwargs...
)
    fnames = Glob.glob(pattern)
    return rotnorth(fnames, offset; kwargs...)
end

function rotnorth(
    fnames::AbstractArray{<:AbstractString},
    offset=[0,0];
    force=false,
    invert=false
)
    outfnames = String[]
    for fname in fnames
        if invert
            outfname = replace(fname,".fits"=>".rotback.fits")
        else
            outfname = replace(fname,".fits"=>".rotnorth.fits")
        end
        if !force && isfile(outfname) && Base.Filesystem.mtime(outfname) > Base.Filesystem.mtime(fname)
            push!(outfnames, outfname)
            continue
        end
        img = load(fname)
        if offset != [0,0]
            img["ROTOFF-X"] = offset[1]
            img["ROTOFF-Y"] = offset[2]
        end
            
        cx, cy = (img["STAR-X"], img["STAR-Y"]) .+ offset

        # TODO: must confirm the rotation direction is correct!
        rot_deg = -img["ANGLE_MEAN"]
        if invert
            rot_deg *= -1
        end
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
        if invert
            push!(applied, History, "$(Date(Dates.now())): Rotated backwards (not North-Up).")
        else
            push!(applied, History, "$(Date(Dates.now())): Rotated North-Up.")
        end
        applied["ORIG.ANGLE_MEAN"] = applied["ANGLE_MEAN"]
        applied["ANGLE_MEAN"] = 0
        AstroImages.writefits(outfname, applied)
        if invert
            println(outfname, "\t rotback ", rot_deg)
        else
            println(outfname, "\t rotnorth ", rot_deg)
        end
        push!(outfnames, outfname)
    end
    return outfnames
end


# Given frames that have been rotated north up, undo the rotation.
function rotrestore(
    pattern::AbstractString,
    offset=[0,0];
    force=false,
)
    fnames = Glob.glob(pattern)
    outfnames = String[]
    for fname in fnames
        outfname = replace(fname,".fits"=>".unrot.fits")
        if !force && isfile(outfname) && Base.Filesystem.mtime(outfname) > Base.Filesystem.mtime(fname)
            push!(outfnames, outfname)
            continue
        end
        img = load(fname)
        if offset != [0,0]
            img["ROTOFF-X"] = offset[1]
            img["ROTOFF-Y"] = offset[2]
        end
            
        cx, cy = (img["STAR-X"], img["STAR-Y"]) .+ offset

        # TODO: must confirm the rotation direction is correct!
        rot_deg = img["ORIG.ANGLE_MEAN"]
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
        push!(applied, History, "$(Date(Dates.now())): Undid north-up rotation.")
        applied["ANGLE_MEAN"] = applied["ORIG.ANGLE_MEAN"]
        AstroImages.writefits(outfname, applied)
        println(outfname, "\t rot undo ", rot_deg)
        push!(outfnames, outfname)
    end
    return outfnames
end
