using Optim
using AstroImages
using ImageFiltering
using CoordinateTransformations, ImageTransformations
using ImageTransformations.Interpolations
using FITSIO
using CairoMakie: Makie
function fluxnorm(pattern)
    fnames = Glob.glob(pattern)
    
    # Load small cutouts around the star
    cutouts = map(fnames) do fname
        img = load(fname)
        x = round(Int, img["STAR-X"])
        y = round(Int, img["STAR-Y"])
        return mapwindow(median, img[x-20:x+20, y-20:y+20], (3,3))
    end

    # Stack and treat this average as the flux target
    target = median(stack(cutouts),dims=3)[:,:]
    display(imview(target))
    target ./= maximum(target)

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
    
    # Now go through images, load 'em up one at a time, and multiply by both
    # factors.
    for (fname,c1,c2) in zip(fnames,c1,c2)
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
    save(replace(pattern, "*"=>"_",".fits"=>"flux.png", ".gz"=>""),figaxpl)
    return c1 .* c2
end