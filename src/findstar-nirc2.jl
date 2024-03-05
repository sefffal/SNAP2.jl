using Dates
using Printf
using Statistics
using FITSIO
using Glob
using Interpolations: Interpolations
using CoordinateTransformations: CoordinateTransformations
using StaticArrays
using OffsetArrays: OffsetArrays
using ImageTransformations: ImageTransformations
using ProgressLogging
using AstroImages
using ImageFiltering

export findstar_nirc2

"""
This function finds the rough location of the star in a sequence of images.
It does this by cross-correlating a template.

The results are returned, and written out to a file.

```
fnames = Glob.glob("caldir.*.2.initial-cal.fits.gz")
load(replace(cal["output_psf_path"], ".fits"=>"-cored.fits"))
findstar_nirc2(fnames, tempalte)
````
"""
function findstar_nirc2(
    fnames, template_fname,
    search_box_px,
    ; verbose=true, force=false, showplots=false)
    verbose && println(SNAP.banner()*"\n\nFIND STAR")

    if fnames isa AbstractString
        fnames = Glob.glob(fnames)
    end


    @info "Performing initial registration to mask star from sky background subtraction:"

    # Get size of initial PSF search region

    # Mostly similar to xcorr_center but doesn't do the translation
    function find_rough_center(
        target, template_img;
        mask=nothing,
        centrefit=true,
        normed=false,
        search_box_px = 160.0,
        showplot=false,
    )
        i1 = round(Int,mean(axes(target,1))-search_box_px)
        i2 = round(Int,mean(axes(target,1))+search_box_px)
        i3 = round(Int,mean(axes(target,2))-search_box_px)
        i4 = round(Int,mean(axes(target,2))+search_box_px)

        

        prepared = target[i1:i2, i3:i4]

        # The cross correlation does not work if there are missing data
        prepared[(!isfinite).(prepared)] .= 0
        prepared = mapwindow(median,prepared,(3,3))
        prepared.=abs.(prepared)
        # display(imview(prepared))


        # The cross correlation does not work if there are missing data
        template = deepcopy(template_img)
        template[(!isfinite).(template)] .= 0
        template.=abs.(template)

        # Optionally apply a mask from a tuple of (min_r, max_r values)
        if !isnothing(mask)
            # Target
            r = sqrt.(prepared.xcoords.^2 .+ prepared.ycoords.^2)
            prepared[(!).(first(mask) .≤ r .≤ last(mask))] .= 0.0

            # Template
            r = sqrt.(template_img.xcoords.^2 .+ template_img.ycoords.^2)
            template[(!).(first(mask) .≤ r .≤ last(mask))] .= 0.0

        end



        z = imfilter(Float64, prepared, OffsetArrays.centered(template),)
        if normed
            q = sqrt.(imfilter(Float64, prepared.^2, OffsetArrays.centered(ones(size(template))),
                           ))
            m = sqrt(sum(template.^2))
            corr = z ./ (m * max.(eps(), q))
        else
            corr = asinh.(z)
        end

        if showplot
            display(imview(corr,clims=extrema,stretch=asinhstretch))
        end

        if centrefit
            # display(imview(corr))
            result = SNAP.centrefit_fixedstd(corr; std=5.0)
            x = result.x
            y = result.y
        else
            # Or just take highest value
            x,y = Tuple(argmax(corr))
        end
        return [x+i1,y+i3,x,y] 
    end

    template = load(template_fname)

    frames = AstroImageMat[]
    fnames_out = String[]
    @info "Loading frames"
    for fname in fnames
        h = FITS(fname,"r") do fits
            read_header(fits[1])
        end
        if force || !haskey(h, "STAR-X")
            push!(frames, load(fname))
            push!(fnames_out, fname)
        end
    end

    # Produce frames
    @info "Finding star"
    tosave = []
    Threads.@threads :dynamic  for (cal,fname,i) in collect(zip(frames,fnames_out,1:length(frames)))
        x,y = find_rough_center(
            cal,
            template,
            centrefit=false,
            normed=true,
            search_box_px=search_box_px,
            showplot=showplots
        )
        push!(cal, History, "$(Date(Dates.now())): Found star.")
        cal["STAR-X"] = x
        cal["STAR-Y"] = y
        push!(tosave, (fname,cal))
        println(fname, "\t($i)\tstar pos: $x\t$y")
        # star_x_y[i] = x,y
    end
    @info "Writing files"
    # wait for all to save or surface errors
    for (fname,cal) in tosave
        AstroImages.writefits(fname,cal)
    end

    # return star_x_y
    return
end