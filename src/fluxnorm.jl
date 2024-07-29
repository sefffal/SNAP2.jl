using Optim
using AstroImages
using ImageFiltering
using CoordinateTransformations, ImageTransformations
using ImageTransformations.Interpolations
using FITSIO
using CairoMakie: Makie


export fluxnorm

function fluxnorm(pattern::AbstractString; kw...)
    fnames = Glob.glob(pattern)
    plotoutfname = replace(pattern, "*"=>"_",".fits"=>"flux.png", ".gz"=>"")
    return fluxnorm(fnames, plotoutfname; kw...)
end
    
function fluxnorm(fnames::AbstractArray{<:AbstractString}, plotoutfname::AbstractString; ignoreinner=5)
        
    # Load small cutouts around the star
    cutouts = map(fnames) do fname
        img = load(fname)
        x = round(Int, img["STAR-X"])
        y = round(Int, img["STAR-Y"])
        cutout = img[x-35:x+35, y-35:y+35] 
        return cutout
    end

    # Stack and treat this average as the flux target
    target = median(stack(cutouts),dims=3)[:,:]
    display(imview(target))
    target ./= maximum(target)

    target[imgsep(AstroImage(target)).<ignoreinner] .= 0


    # Find the factor we need to multiply each frame by to reach the average
    c1 = map(cutouts) do img
        vec(img) \ vec(target)
    end
    cutouts = broadcast.(*, cutouts, c1) #map((img,c)->img .* c, 

    # Stack and treat this average as the flux target (round 2)
    target = median(stack(cutouts),dims=3)[:,:]

    # Find the factor we need to multiply each frame by to reach the average
    c2 = map(cutouts) do img
        vec(img) \ vec(target)
    end

    if any(!isfinite, c1) || any(!isfinite, c1)
        c1 .= c2 .= 1/length(c1)
    end
    
    # Now go through images, load 'em up one at a time, and multiply by both
    # factors.
    for (fname,c1,c2) in zip(fnames,c1,c2)
        if !isfinite(c1*c2)
            @warn "non finite coefficients" c1 c2 fname
            c1 = c2 = 0
        end
        img = load(fname).*c1.*c2
        push!(img, History, "$(Date(Dates.now())): Flux normalized to median of sequence.")
        img["FLUXNORM"] = c1 * c2 
        img["FLUXNORM",Comment] = "image multiplied by this flux normalization factor"

        outname = replace(fname, ".fits"=>".fluxnorm.fits",)
        AstroImages.writefits(outname, img)
        println(outname, "\t flux ", c1*c2)
    end

    figaxpl = Makie.lines(
        1 ./ (c1 .* c2),
        axis = (;
            xlabel = "frame #",
            ylabel = "peak flux [counts]",
        )
    )
    save(plotoutfname,figaxpl)
    return c1 .* c2
end