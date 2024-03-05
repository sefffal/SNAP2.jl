using Glob
using ImageTransformations
using CoordinateTransformations: CoordinateTransformations
using Rotations: Rotations
using LinearAlgebra


export rotnorthrefs

"""
Given a pattern matching multiple calibrated, registered files, 
rotate all frames North-up, and prepare a folder of all reference images.
Reference images will be rotated and/or scaled such their speckles are still
aligned with the target image.
"""
function rotnorthrefs(pattern::AbstractString, refpattern::AbstractString=pattern; force=false)
    fnames = Glob.glob(pattern)
    fnames_refs = Glob.glob(refpattern)

    for fname in fnames
        # Rotate north or get filename if already done.
        # rotnorthfname = only(rotnorth(fname))
        fs = rotnorth(fname; force)
        rotnorthfname = only(fs)
        targ = load(rotnorthfname)
        refdir = replace(rotnorthfname, ".fits.gz"=>".refs")
        if !isdir(refdir)
            mkdir(refdir)
        end
        for fname_ref in fnames_refs
            outfname = joinpath(refdir, splitpath(fname_ref)[end])
            if !force && isfile(outfname) && Base.Filesystem.mtime(outfname) > Base.Filesystem.mtime(fname)
                continue
            end

            ref = load(fname_ref)
            cx, cy = ref["STAR-X"], ref["STAR-Y"]

            #  Undo our rotation, then do their rotation so that we match
            rot_deg = -targ["ORIG.ANGLE_MEAN"]
            rot_rad = deg2rad(rot_deg)
            tfrm = CoordinateTransformations.LinearMap(
                Rotations.RotMatrix(rot_rad)
            )

            # Create a coordinate-transform that maps the data into
            # speckle-aligned coordinates (pupil tracking in ADI)
            if ref["WAVE"] != targ["WAVE"]
                tfrm = LinearMap(UniformScaling(targ["WAVE"] / ref["WAVE"] )) ∘ tfrm
            end

            tfrm = CoordinateTransformations.recenter(tfrm, (cx,cy))
            applied_dat = collect(ImageTransformations.warp(
                collect(ref),
                tfrm,
                axes(ref),
                # Broken for some reason:...
                # Interpolations.Quadratic(Interpolations.Flat(Interpolations.OnGrid())),
                NaN,
            ))
            applied = copyheader(ref, applied_dat)
            push!(applied, History, "$(Date(Dates.now())): Rotated to match target img North-up angle.")
            applied["ANGLE_MEAN"] = ref["ANGLE_MEAN"] - targ["ORIG.ANGLE_MEAN"]
            AstroImages.writefits(outfname, applied)
            println(outfname, "\t rot ref ", rot_deg)
        end
    end
end