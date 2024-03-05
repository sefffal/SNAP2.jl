
function header_normalize!(
    img::AstroImage,
)

    # TODO: update

    angles = Float64[]
    if haskey(img, "PA_DETLEN")
        for i in 1:img["PA_DETLEN"]
            push!(angles, img[(@sprintf "PA_DET%0d" i)])
        end
    end
    if length(angles) < 1 && haskey(img, "PA") && img["PA"] != 0.7
        push!(angles, img["PA"])
    end
    # Override to ensure this is correct.
    # Some codes (including old versions) do not handle averaging angles
    # correctly.
    if length(angles) > 0
        img["ANGLE_MEAN"] = meandegrees(angles)
    end

    # Sometimes keys are wrapped in strings, and sometimes not.
    # This will work with either.
    tryparseconvert(str::String) = parse(Float64, str)
    tryparseconvert(flt::Float64) = flt
    tryparseconvert(num::Number) = convert(Float64, num)

    # Handle data scaling
    if haskey(img, "BSCALE")
        img .*= tryparseconvert(img["BSCALE"])
    end
    if haskey(img, "BZERO")
        img .+= tryparseconvert(img["BZERO"])
    end

    # if haskey(img, "TOTEXP")
    #     img["TOTEXP"] = img["TOTEXP"] * 60 * 60
    #     set_comment!(img, "TOTEXP", "Exposure in seconds")
    # elseif haskey(img, "ITIME") && haskey(img, "COADDS")
    #     img["TOTEXP"] =
    #         Float64(img["ITIME"]) * Float64(img["COADDS"])
    #     set_comment!(img, "TOTEXP", "Exposure in seconds")
    # end

    # Set up the axes correctly using header, with fallbacks
    xlen = size(img, 1)
    ylen = size(img, 2)
    if haskey(img, "PA_DETLEN")
        if haskey(img, "CDELT2")
            xscale = tryparseconvert(img["CDELT1"]) * -1e3 * 3600
            yscale = tryparseconvert(img["CDELT2"]) * 1e3 * 3600
        elseif haskey(img, "PIXSCAL1") && haskey(img, "PIXSCAL2")
            xscale = tryparseconvert(img["PIXSCAL1"]) * -1e3 * 3600
            yscale = tryparseconvert(img["PIXSCAL2"]) * 1e3 * 3600
        elseif haskey(img, "PIXSCALE")
            xscale = yscale = tryparseconvert(img["PIXSCALE"]) * 1e3 * 3600
            xscale = -1 .* xscale
        else
            xscale = yscale = 1
        end
    else
        if haskey(img, "CDELT1") && haskey(img, "CDELT2")
            # xscale = img["CDELT1"] * -1e3 * 3600
            xscale = tryparseconvert(img["CDELT1"]) * 1e3
            # yscale = img["CDELT2"] * 1e3 * 3600
            yscale = tryparseconvert(img["CDELT2"]) * 1e3
        elseif haskey(img, "PIXSCAL1") && haskey(img, "PIXSCAL2")
            # xscale = img["PIXSCAL1"] * -1e3 * 3600
            xscale = tryparseconvert(img["PIXSCAL1"]) * 1e3
            # yscale = img["PIXSCAL2"] * 1e3 * 3600
            yscale = tryparseconvert(img["PIXSCAL2"]) * 1e3
        elseif haskey(img, "PIXSCALE")
            # xscale = yscale = img["PIXSCALE"] * 1e3 * 3600
            xscale = yscale = tryparseconvert(img["PIXSCALE"]) * 1e3
            xscale = -1 .* xscale
        else
            xscale = yscale = 1
        end
    end

    if haskey(img, "CRPIX1") && haskey(img, "CRPIX2")
        xzero = tryparseconvert(img["CRPIX1"])
        yzero = tryparseconvert(img["CRPIX2"])
    else
        # Fall back to assuming the centre pixel
        xzero = xlen / 2 + 0.5
        yzero = ylen / 2 + 0.5
    end

    if haskey(img, "CRVAL1") && haskey(img, "CRVAL2")
        RA = haskey(img, "RA") ? tryparseconvert(img["RA"]) : 0.0
        DEC = haskey(img, "DEC") ? tryparseconvert(img["DEC"]) : 0.0
        xshift = (tryparseconvert(img["CRVAL1"]) - RA) * 1e3 * 60 * 60 # Convert from degrees to mas
        yshift = (tryparseconvert(img["CRVAL2"]) - DEC) * 1e3 * 60 * 60 # Convert from degrees to mas
    else
        # Fall back to assuming the centre pixel
        xshift = 0.0
        yshift = 0.0
    end

    # xrange = ((1.0:xlen) .- xzero) .* xscale .+ xshift
    # yrange = ((1.0:ylen) .- yzero) .* yscale .+ yshift

    img["CRPIX1"] = xzero
    img["CRPIX2"] = yzero
    img["CRVAL1"] = xshift
    img["CRVAL2"] = yshift
    img["CDELT1"] = xscale
    img["CDELT2"] = yscale

    # TODO: see if necessary
    # for I in 1:length(img)
    #     #@warn "Replacing NaN header value with nothing" key=img.keys[I] maxlog=10
    #     # Commented out this warning. It causes a lot of warning spam in the logs.
    #     if typeof(img.values[I])<:AbstractFloat && isnan(img.values[I])
    #         img.values[I] = nothing
    #     end
    # end

    # TODO: see if necessary?
    # # Handle position angle
    # if !haskey(img, "PA")
    #     mean_PA = SNAP.meandegrees(frame.angles)

    #     img["PA"] = rad2deg(rem2pi(deg2rad(mean_PA), RoundDown))
    # end
    # img["PA_DETLEN"] = length(frame.angles)
    # set_comment!(img, "PA_DETLEN", "Number of detailed PA values to follow")
    # for (i, a) in pairs(frame.angles)
    #     key = (@sprintf "PA_DET%0d" i)
    #     img[key] = rad2deg(rem2pi(deg2rad(a), RoundDown))
    #     set_comment!(img, key, "Telescope position angle during second $i")
    # end

    # Additional telescope specific processing
    if header_getinst(img) == :NIRC2
        header_normalize_nirc2!(img)
    end

    return img
end

function header_getinst(img::AstroImage)
    if haskey(img, "TELESCOP") &&
        startswith(img["TELESCOP"], "Keck") &&
        haskey(img, "CURRINST") &&
        img["CURRINST"] == "NIRC2"
        return :NIRC2
    else
        @warn "instrument not detected or implemented"
    end
end


function header_setpipelinehash!(img::AstroImage)
    # Try to put the git commit hash of the pipeline into the img of each file
    # to help with reproducibility
    try
        gitrepo = joinpath(dirname(pathof(SNAP)), "..")
        hash = string(LibGit2.head(gitrepo))
        if LibGit2.isdirty(LibGit2.GitRepo(gitrepo))
            hash *= "-dirty"
        end
        img["PIPELINE"] = hash
        img["PIPELINE",Comment] = "Git commit hash of the code used to produce this image. -dirty if there were uncommitted working copy changes."
        @info "Saving pipeline commit hash in img" hash maxlog=1
    catch err
        @warn "Could not store pipeline version/commit hash in header" exception=err
    end

    return img
end


function header_set_pa_angles!(img::AstroImage, angles::AbstractArray{<:Number})
    img["PA_DETLEN"] = length(angles)
    img["PA_DETLEN",Comment] = "Number of PA values to follow"
    for (i,a) = pairs(frame.angles)
        key = @sprintf("PA_DET%0d", i)
        img[key] = a
        img[key,Comment] = "Telescope position angle during second $i"
    end
    return img
end
function header_get_pa_angles(img::AstroImage)
    angles = Float64[]
    if haskey(img, "PA_DETLEN")
        for i=1:parse(Int,headers["PA_DETLEN"])
            push!(angles, headers[(@sprintf "PA_DET%0d" i)])
        end
    end
    return angles
    return img
end


# Correctly average angle measurements
meandegrees(degrees) = if isempty(degrees)
    zero(eltype(degrees))
else
    rad2deg(atan(sum(sind.(degrees)), sum(cosd.(degrees))))
end
meanrad(rads) = if isempty(rads)
    zero(eltype(rads))
else
    atan(sum(sin.(rads)), sum(cos.(rads)))
end

include("header-normalize-nirc2.jl")