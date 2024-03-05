#=
This script prepares frames for subtraction.

One script rotates images North-up.

Another script rotates and/or scales reference frames to match a given target frame.

=#

export inject

using LinearAlgebra
function inject(
    pattern::AbstractString;
    pa_deg,
    sep_px,
    force=false,
    psf_fwhm,
    psf_amp=1.0,
    psf_ratio=0.0,
    invert=false
)
    fnames = Glob.glob(pattern)
    outfnames = String[]
    i = 0
    for fname in fnames
        i+= 1
        if invert
            outfname = replace(fname,".fits"=>".injectback.fits")
        else
            outfname = replace(fname,".fits"=>".inject.fits")
        end
        if !force && isfile(outfname) && Base.Filesystem.mtime(outfname) > Base.Filesystem.mtime(fname)
            push!(outfnames, outfname)
            continue
        end
        img = load(fname)
            
        # cx, cy = (img["STAR-X"], img["STAR-Y"])

        # TODO: must confirm the rotation direction is correct!
        rot_deg = -img["ANGLE_MEAN"]
        if invert
            rot_deg *= -1
        end
        rot_rad = deg2rad(rot_deg + pa_deg)

        x = sep_px * cos(rot_rad)
        y = sep_px * sin(rot_rad)
        psf_model = PSFModels.airydisk(;
            amp   = psf_amp,
            fwhm  = psf_fwhm,
            ratio = psf_ratio,
            x,y
        )
        seps = imgsep(img)
        thetas = imgang(img)
        applied =  img .+ psf_model.(
            seps .* cos.(thetas),
            seps .* sin.(thetas)
        )

        if invert
            push!(applied, History, "$(Date(Dates.now())): Injected planet backwards")
        else
            push!(applied, History, "$(Date(Dates.now())): Injected planet.")
        end
        applied["INJ_PA"] = pa_deg
        applied["INJ_PA",Comment] = "PA [deg] of injected planet"
        applied["INJ_SEP"] = sep_px
        applied["INJ_SEP",Comment] = "Sep [px] of injected planet"
        AstroImages.writefits(outfname, applied)
        if invert
            println(outfname, " $(i) \t inject back ", rad2deg(rot_rad))
        else
            println(outfname, " $(i) \t inject ", rad2deg(rot_rad))
        end
        push!(outfnames, outfname)
    end
    return outfnames
end
