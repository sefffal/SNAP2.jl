export matchedfilt
# Given a file that has been PSF subtracted, whose second FITS HDU is a 
# cube of PSF models associated with each region
function matchedfilt(
    pattern::AbstractString,
    regions_S::AbstractVector;
    force=false
)
    fnames = Glob.glob(pattern)
    for fname in fnames
        outfname = replace(fname,".fits"=>".mf.fits")
        if !force && isfile(outfname) && Base.Filesystem.mtime(outfname) > Base.Filesystem.mtime(fname)
            continue
        end
        img = load(fname,1)
        models = load(fname,2)
        if size(models,3) != length(regions_S)
            error(
                "Number of regions provided do not match the number of output PSF models in the provided image. "*
                "Maybe the regions are not the same as the ones used in the reduction?"
            )
        end
        out = deepcopy(img)
        fill!(out, NaN)

        # Get contrast curve
        sep, cont = only(contrast(fname; force))
        contrast_interpolator = Interpolations.linear_interpolation(
            sep, cont,
            extrapolation_bc=Line()
        )

        matchedfilt!(out, img, regions_S, models, contrast_interpolator)

        push!(out, History, "$(Date(Dates.now())): Applied matched-filter.")
        AstroImages.writefits(outfname, out)
        println(outfname, "\t mf ")
    end
end


function matchedfilt!(
    out::AstroImage,
    img::AstroImage,
    regions_S::AbstractVector,
    models::AbstractArray,
    contrast_interpolator
)
    rs = imgsep(img)
    θs = imgang(img)
    xc = mean(axes(img,1))
    yc = mean(axes(img,2))

    model_interpolators = map(axes(models,3)) do model_i
        m = models[:,:,model_i]
        m[.!isfinite.(m)] .= 0
        Interpolations.linear_interpolation(
            (axes(models,1), axes(models,2)),
            m,
            extrapolation_bc=0
        )
    end

    # Determine centre location of the PSF model in each region
    # TODO: better to store this as a header or something during subtraction and then just read it back here.
    region_θs = Float64[]
    region_rs = Float64[]
    for region_S in regions_S
        inject_planet_sep_ave = mean(rs[region_S.>0])
        inject_planet_pa_ave = deg2rad.(SNAP.meandegrees(vec(rad2deg.(θs[region_S.>0]))))
        push!(region_θs, inject_planet_pa_ave)
        push!(region_rs, inject_planet_sep_ave)
    end

    # We now do the convolution. Our kernel varies as a function of position, so we have to do 
    # this ourselves naively looping over all pixels.

    # Go through each output pixel of the image (pixel != NaN)
    # I_out is the location of the pixel we will output in this iteration of the loop.
    Threads.@threads  for I_out in findall(isfinite, img) # Returns CartesianIndices

        # Total value of kernel and size of kernel for normalization.
        ktot = 0.0
        kernsize = 0
        out[I_out] = 0

        # Consider that for this output pixel location, we have a PSF model centred here.
        # We then go over all pixels and add up how much they agree with that PSF model.
        output_pixel_r = rs[I_out]
        output_pixel_θ = θs[I_out]

        # output_pixel_r = 40.6
        # output_pixel_θ = deg2rad(42)



        # Loop through *all* finite pixels in the input image to consider their contribution
        # to a hypothetical planet at this output location.
        for I_in in findall(isfinite, img)
            # Location of this input pixel
            input_pixel_r = rs[I_in]
            input_pixel_θ = θs[I_in]

            # Different input pixels also have different forward models (each subtraction region
            # has its own model). We loop through the different pairs of (region, model) to select the
            # correct one.
            # Assumption: the subtraction regions do not overlap.
            the_i_model = nothing
            for i_model in axes(models,3)
                if regions_S[i_model][I_in] > 0
                    the_i_model = i_model
                    break
                end
            end
            # If no input model for this pixel, warn and continue.
            # Something probably went wrong since we intend to always output a model
            if isnothing(the_i_model) 
                @warn "no matching PSF model for this input pixel (maxlog=1)" I_in maxlog=1
                continue
            end

            # This is the location of the centre of the planet model for this region.
            model_cent_r = region_rs[the_i_model]
            model_cent_θ = region_θs[the_i_model]

            # Now we have to interpolate into the correct location of the model.
            # Consider 1st how we have to translate the model to our desired centration on the output 
            # pixel:
            # If we are looking at output_pixel_r + Δ, we need to instead index into `model_cent_r + Δ`
            # where Δ is our displacement between input and output pixels (positive as input_pixel_r > output_pixel_r).
            Δ_r = output_pixel_r - input_pixel_r
            Δ_θ = output_pixel_θ - input_pixel_θ

            x = xc + (model_cent_r + Δ_r) * cos(model_cent_θ + Δ_θ)
            y = yc + (model_cent_r + Δ_r) * sin(model_cent_θ + Δ_θ)



            # Note: it might be more efficient to invert these loops?
            model_at_I_in = model_interpolators[the_i_model](x,y)
            if model_at_I_in == 0.0 # Exactly zero
                continue
            end
            noise = contrast_interpolator(input_pixel_r)
            if !isfinite(noise)
                @warn "non finite contrast where finite was expected" maxlog=1
            end
            if noise <= 0
                @warn "zero valued or negative contrast will lead to bad matched filter output" maxlog=1
            end
            model_at_I_in /= noise^2
            out[I_out] += img[I_in]*model_at_I_in
            ktot += model_at_I_in^2
            kernsize += 1
        end
        out[I_out] /= sqrt(ktot/kernsize)
        if !isfinite(out[I_out])
            @warn "non finite matched filter output where finite was expected" maxlog=1
        end
    end

    return out
end


    # # for now don't consider scaling due to SDI, only consider rotation from ADI.
    # x = inject_planet_sep_ave * cos(inject_planet_pa_ave)
    # y = inject_planet_sep_ave * sin(inject_planet_pa_ave)

    # Matched filter
    # Correlate each pixel with the template PSF at that location.
    
    # How would this work? We don't have a 2D image but a 1D vector of values
    # and their positions.
    
    # For each output pixel
        # Call this I
        # Loop through each spatial location
            # Call this J
            # k := Calculate FM of PSF centred at I at location J
            # add k to kernel for normalization
            # multiply k by output pixel value
    
    # for I in findall(region_S)
#         data_fmmfout[I] = 0
#         pos_I = SVector(Sx[I], Sy[I])
#         ktot = 0.0
#         kernsize = 0
#         for J in eachindex(data_out, data_fmmfout, Sx, Sy)
#             pos_J = SVector(Sx[J], Sy[J])
#             sep = sqrt((Sx[I]-Sx[J])^2+(Sy[I]-Sy[J])^2)
#             # TODO: this is a hardcoded upper limit on convolution kernel size  for performance
#             if sep > 200 
#                 continue
#             end
#             r = sqrt(Sx[I]^2+Sy[I]^2)
#             noise = contrast_curve(r)
#             if !isfinite(noise)
#                 continue
#             end
#             k = 0.0
#             for (r,c) in zip(references, coefficients)
#                 k += c*Companion.fm_phot(r, pos_J, pos_I)
#             end
#             k /= noise^2
#             data_fmmfout[I] += data_out[J]*k
#             ktot += k^2
#             kernsize += 1
#         end
#         data_fmmfout[I] /= sqrt(ktot/kernsize)
#     end


# end