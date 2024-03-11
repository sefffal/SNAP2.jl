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

        @info "filtering "
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
    # TODO: do we need to do the filtering across multiple images at once? Or do we just add them together afterwards?

    # TODO: need to get the contrast curve too for noise-weighting.

    # Determine centre location of the PSF model.
    # TODO: better to store this as a header or something during subtraction and then just read it back here.


    rs = imgsep(img)
    θs = imgang(img)
    xc = mean(axes(img,1))
    yc = mean(axes(img,2))


    # @warn "debug workaround bad pixels"
    # m = .!isfinite.(img)
    # img[begin:651,begin:651] .= 0
    # img[m] .= NaN

    # Location of model centre for each region
    region_θs = Float64[]
    region_rs = Float64[]
    for region_S in regions_S
        inject_planet_sep_ave = mean(rs[region_S.>0])
        inject_planet_pa_ave = deg2rad.(SNAP.meandegrees(vec(rad2deg.(θs[region_S.>0]))))
        push!(region_θs, inject_planet_pa_ave)
        push!(region_rs, inject_planet_sep_ave)
    end

    # Go through each output pixel of the image (pixel != NaN)
    Threads.@threads  for I_out in findall(isfinite, img) # Returns CartesianIndices
        # Total value of kernel, for normalization.
        ktot = 0.0
        kernsize = 0
        out[I_out] = 0
        # Imagine, for this output pixel location, that we have a PSF model centred here.
        # We then go over all pixels and add up how much they agree with that PSF model.

        # The tricky parts are 1) handling rotation and separation of the model vs this location.
        # [we could just store coefficients and re-calculate here?] and 2) handling how the regions
        # split things up, as they all have different models.
        
        # This gives the coordinates into the models for a model centred on this output pixel,
        # accounting for the model's original position.
        r_model_wantcent = rs[I_out] #region_rs[i_model] + rs[I_out]
        θ_model_wantcent = θs[I_out] #region_θs[i_model] + θs[I_out]


        # For handling region seams each with different models: we have to select the right model
        # based on the input pixel we are currently examining.
        for I_in in findall(isfinite, img)
            this_r = rs[I_in]
            this_θ = θs[I_in]
            # Select correct model.
            the_i_model = nothing
            # At this input pixel, there is exactly one matching model. Find it.
            for i_model in axes(models,3)
                if regions_S[i_model][I_in] > 0
                    the_i_model = i_model
                    break
                end
            end
            if isnothing(the_i_model) 
                continue
            end
            r_model_wascent = region_rs[the_i_model]
            θ_model_wascent = region_θs[the_i_model]

            # x = (r_model_wantcent)*cos(θ_model_cent+δθ)
            # y = (r_model_cent+δr)*sin(θ_model_cent+δθ)
            # @show r_model_cent rs[I_out] rs[I_in] region_rs[the_i_model] 


            # Need to look into model at (was cent) + our offset, which is (this_ - wantcent)
            # Add our displacement vs where we centred the model to the original location
            x = round(Int, xc + (r_model_wascent + this_r - r_model_wantcent)*cos(θ_model_wascent + this_θ - θ_model_wantcent))
            y = round(Int, yc + (r_model_wascent + this_r - r_model_wantcent)*sin(θ_model_wascent + this_θ - θ_model_wantcent))
            # x = round(Int, xc + (r_model_wascent - this_r + r_model_wantcent)*cos(θ_model_wascent - this_θ + θ_model_wantcent))
            # y = round(Int, yc + (r_model_wascent - this_r + r_model_wantcent)*sin(θ_model_wascent - this_θ + θ_model_wantcent))

            # This gives x and y, but not the actual indices. For that we'll need to add the centre index of the image.

            
            # Note: it might be more efficient to invert these loops?
            # if regions_S[the_i_model][x,y]>0 && 
            model_at_I_in = nothing
            if isfinite(models[x,y,the_i_model])
                # TODO: we will have to shift the model around here.
                # model_at_I_in = models[I_in,the_i_model]
                model_at_I_in = models[x,y,the_i_model] # TODO: interpolate
            end
            if isnothing(model_at_I_in) 
                continue
            end
            noise = contrast_interpolator(this_r)
            if !isfinite(noise)
                @warn "non finite contrast where finite was expected" maxlog=1
            end
            if noise <= 0
                @warn "zero valued or negative contrast will lead to bad matched filter output" maxlog=1
            end
            # TODO: need noise here form contrast curve
            model_at_I_in /= noise^2
            # model_at_I_in /= noise
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