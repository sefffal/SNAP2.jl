using Random: Random
using Enzyme
using Optim: Optim

# TODO:
# At snropt_all level:
#     * Add parameters for chunk size and step.

# Models
# Output PSF models too. Store these as a cube in a second extension to the same output image.
# We really want to compute the model on a larger area than just the subtraction region.
# We need a third region that is the union of the subtraction, optimization, and buffers.
# Should check if it's fast to just calculate over the full frame.

# FMMF
# Post processing step 
# Can I do this right on the subtracted image, despite the discontinuities?
# I think so!
# just need to be picking from the right model region.
# I can also rotate the model as needed, instead of assuming it doesn't change based on PA. The negative wings definitely do!.
# Estimating the noise: need to do a noise-weighted matched filter.
# But how do we do esimate the noise? I guess we can do a contrast curve
#     from the subtracted image and assume it's consistent between annular sectors.
# Could also do this on the full stack if we used a mean, but not easily if we used a median.

export snropt_region!, snropt_frame, snropt_all, snropt_multitarg

function snropt_all(fnames_pattern, regions_S, regions_O, regions_M; force=false, kwargs...)
    fnames = glob(fnames_pattern)
    imgs = load.(fnames)
    refcube=stack(imgs)
    i = 0
    for (fname, targ) in zip(fnames, imgs)
        i+=1
        if haskey(targ, "ISRDIREF") && targ["ISRDIREF"] > 0
            println(fname, "\t($i) skipping reference")  
            continue
        end
        outfname = replace(fname, ".fits"=>".snropt.fits")
        if !force && isfile(outfname) &&  Base.Filesystem.mtime(outfname) > Base.Filesystem.mtime(fname)
            out = load(outfname,1)
        else
            out = deepcopy(targ)
            fill!(out, NaN)
        end
        fm = similar(out, size(out)..., length(regions_S))
        fill!(fm, NaN)

        # for (reg_S, reg_O) in zip(regions_S, regions_O)
        # Threads.@threads :dynamic
        i = 0 
        for (reg_S, reg_O, reg_M) in collect(zip(regions_S, regions_O, regions_M))
            i += 1
            # TODO: need a way to know where each model is centred. Might not be peak! Or is?
            snropt_region!(out, view(fm, :,:,i), targ, imgs, reg_S.>0, reg_O.>0, reg_M.>0, refcube; kwargs...)
        end
        AstroImages.writefits(outfname, out, fm)      
        println(outfname, "\t($i)")  
    end
end

"""
This function sets up and performs a multi-target SNR optimization.

"""
function snropt_multitarg(
    fnames_pattern,
    regions_S,
    regions_O,
    regions_M;
    force=false,
    target_group=10,
    target_step=10,
    target_start=1,
    target_stop=length(glob(fnames_pattern)),
    save_refs=false,
    kwargs...
)
    fnames = glob(fnames_pattern)

    for target_i_start in target_start:target_step:target_stop
        target_i_range = target_i_start:1:min(target_i_start + target_group - 1, lastindex(fnames))

        # Start by preparing the reference images for each of the targets
        all_refs = AstroImageMat[]
        for target_i in target_i_range
            these_refs = rotnorthrefs(fnames[target_i], fnames_pattern; force, save=save_refs)
            append!(all_refs, these_refs)
        end

        targ_fname = replace(fnames[target_i_start], ".fits"=>".rotnorth.fits")

        # Do the subtraction for these targets
        snropt_frame(
            # Target image (we actually have multiple targets, but this sets eg the rotation angle for modelling)
            targ_fname,
            # Reference images rotated to match the north-up target images
            all_refs,
            regions_S,
            regions_O,
            regions_M;
            force,
            kwargs...
        )

    end
end


#     i = 0
#     for (fname, targ) in zip(fnames, imgs)
#         i+=1
#         outfname = replace(fname, ".fits"=>".snropt.fits")
#         if !force && isfile(outfname) &&  Base.Filesystem.mtime(outfname) > Base.Filesystem.mtime(fname)
#             out = load(outfname)
#         else
#             out = deepcopy(targ)
#             fill!(out, NaN)
#         end

#         # for (reg_S, reg_O) in zip(regions_S, regions_O)
#         # Threads.@threads :dynamic 
#         for (reg_S, reg_O) in collect(zip(regions_S, regions_O))
#             snropt_region!(out, targ, imgs, reg_S.>0, reg_O.>0, refcube; kwargs...)
#         end
#         AstroImages.writefits(outfname, out)      
#         println(outfname, "\t($i)")  
#     end
# end

# # This is an implementation that works on multiple cores
# # To use it, you must first start julia with multiple workers (`julia -p 8`) and
# # then run `using Distributed`, `@everywhere using SNAP`.
# function snropt_all_distributed(fnames_pattern, regions_S, regions_O; force=false, kwargs...)
#     fnames = glob(fnames_pattern)
#     @everywhere imgs = load.(fnames)
#     @everywhere refcube=stack(imgs)
#     pmap(eachindex(fnames)) do i
#         fname = fnames[i]
#         outfname = replace(fname, ".fits"=>".snropt.fits")
#         if !force && isfile(outfname) &&  Base.Filesystem.mtime(outfname) > Base.Filesystem.mtime(fname)
#             out = load(outfname)
#         else
#             out = deepcopy(targ)
#             fill!(out, NaN)
#         end
#         for (reg_S, reg_O) in collect(zip(regions_S, regions_O))
#             snropt_region!(out, targ, imgs, reg_S.>0, reg_O.>0, refcube; kwargs...)
#         end
#         AstroImages.writefits(outfname, out)      
#         println(outfname, "\t($i)")
#     end
# end

"""
Given a glob pattern string, run `Glob.glob` and return a vector of filenames.
Given a vector of globs, run `Glob.glob` on each and return a concatenated
vector of filenames.
"""
globvec(pattern::AbstractString) = Glob.glob(pattern)
globvec(patterns::AbstractVector{<:AbstractString}) = reduce(vcat, Glob.glob.(patterns))
function snropt_frame(fname, refnames_pattern, regions_S, regions_O, regions_M;force=false, kwargs...)
    target = load(fname)
    outfname = replace(fname, ".fits"=>".snropt.fits")
    @show force
    if !force && isfile(outfname) && Base.Filesystem.mtime(outfname) > Base.Filesystem.mtime(fname)
        out = load(outfname)
    else
        out = deepcopy(target)
        fill!(out, NaN)
    end

    if refnames_pattern isa AbstractString
        refnames = globvec(refnames_pattern)
        @info "Loading images"
        @progress refs = [
            load(refname)
            for refname in refnames
        ]
    elseif refnames_pattern isa AbstractArray
        refs = refnames_pattern
    end
    # For processing, we want the data in a single large cube
    # and the headers. Once we grab the headers and data julia can 
    # GC the images themselves.
    hdrs = header.(refs)
    @info "Stacking images into contiguous memory"
    refcube = stack(refs)
    fm = similar(out, size(out)..., length(regions_S))
    fill!(fm, NaN)
    i = 0
    for (reg_S, reg_O, reg_M) in zip(regions_S, regions_O, regions_M)
        i+= 1
        println("region #$i")
        snropt_region!(
            view(out, :,:),
            view(fm, :,:, i), 
            target,
            hdrs,
            reg_S.>0,
            reg_O.>0,
            reg_M.>0,
            refcube;
            kwargs...
        )
    end
    AstroImages.writefits(outfname, out, fm)
    println(outfname, "\t (snropt)")
    return out
end



function snropt_region!(
    out::AstroImage, # Subtracted output (only inside regionS)
    fm::AbstractArray,  # Forward modell output (potentially full image)
    target::AstroImage,
    refs::AbstractVector,
    region_S::AbstractMatrix{Bool},
    region_O::AbstractMatrix{Bool},
    region_M::AbstractMatrix{Bool}=region_S,
    refcube=stack(refs);
    N_SVD=0,
    psf_fwhm,
    psf_ratio=0.0,
)
    if count(region_S) == 0 || count(region_O) == 0
        @warn "cannot perform PSF subtraction with an empty region " count(region_S) count(region_O)
    end

    rs = imgsep(target)
    θs = imgang(target)


    refs_are_rdi = map(refs) do ref
        if haskey(ref, "ISRDIREF") && ref["ISRDIREF"] > 0
            return 1
        end
        return 0
    end

    # Given a list of frames, calculate the photometry of the planet at the mean location
    # of the target frame.
    # We can use that for throughput normalization, for tuning hyper-parameters, and for 
    # SNR optimization.

    # We simulate the planet at the average PA and average separation of the target subtraction
    # region.
    inject_planet_sep_ave = mean(rs[region_S])
    inject_planet_pa_ave = deg2rad.(SNAP.meandegrees(vec(rad2deg.(θs[region_S]))))
    # We now figure out how much the planet has rotated in each reference frame.
    # TODO: confirm the sign is correct in second part of expression (no effect if PSF is symmetric though.)
    ref_planet_pa = inject_planet_pa_ave .+ deg2rad.(getindex.(refs, "ANGLE_MEAN").- target["ANGLE_MEAN"] )
    ref_planet_pa .= rem2pi.(ref_planet_pa, RoundDown)

    # for now don't consider scaling due to SDI, only consider rotation from ADI.
    x = inject_planet_sep_ave * cos(inject_planet_pa_ave)
    y = inject_planet_sep_ave * sin(inject_planet_pa_ave)
    psf_model = PSFModels.airydisk(;
        amp   = 1.0,
        fwhm  = psf_fwhm,
        ratio = psf_ratio,
        x=x+eps(x),y=y+eps(y)
        # Eps to work around this bug: https://github.com/JuliaAstro/PSFModels.jl/issues/14
    )

    # Calculate just the point photometry at the peak of the target
    phot_of_target_at_inject_location = psf_model(x,y)
    phot_of_ref_at_target_inject_location = psf_model.(
        inject_planet_sep_ave .* cos.(ref_planet_pa),
        inject_planet_sep_ave .* sin.(ref_planet_pa),
    )
    phot_of_ref_at_target_inject_location .*= (1 .- refs_are_rdi)

    refs_O = zeros(eltype(refcube), count(region_O), length(refs))
    refs_S = zeros(eltype(refcube), count(region_S), length(refs))
    for i in eachindex(refs)
        refs_O[:,i] .= refcube[region_O,i]
        refs_S[:,i] .= refcube[region_S,i]
    end
    valid = vec(all(isfinite,refs_O,dims=2))

    # Reject pixels that are strongly deviated from the clipped standard deviation.
    clip_std = std(StatsBase.trim(vec(view(refs_O,valid,:)); prop=0.10)) # Find std, while ignoring 10% most deviated pixels
    valid .&= all(
        (-4clip_std .< view(refs_O,valid,:) .< 4clip_std),
        dims=2
    ) 
    # Don't do this clipping if we end up with nothing left, or take too much
    if count(valid) == 0|| mean(valid) < 0.25
        valid .= true
    end
    refs_O = refs_O[valid,:]

    # Identify the pixels that are safe to use for evaluating the SNR of the subtraction
    # Avoid looking at bad pixels
    # Ignore the bottom and top 0.5% of pixels
    # clip_ex = extrema(StatsBase.trim(vec(refs_S); prop=0.01)) # Find std, while ignoring 5% most deviated pixels
    # region_S_valid = all(
    #     (clip_ex[1] .< refs_S .< clip_ex[2]),
    #     dims=2
    # )[:]
    # clip_std = std(StatsBase.trim(vec(refs_S); prop=0.10)) # Find std, while ignoring 10% most deviated pixels
    # region_S_valid = all(
    #     (-4clip_std .< refs_S .< 4clip_std),
    #     dims=2
    # )[:]
    # # Don't do this clipping if we end up with nothing left, or take too much
    # if count(region_S_valid) == 0 || mean(region_S_valid) < 0.25
    #     region_S_valid .= true
    # end


    # Calculate models consisting only of the contamination-free synthetic planet PSF
    refs_M = zeros(eltype(refcube), count(region_M), length(refs))
    for i in eachindex(refs)
        x = inject_planet_sep_ave * cos(ref_planet_pa[i])
        y = inject_planet_sep_ave * sin(ref_planet_pa[i])
        psf_model = PSFModels.airydisk(;
            amp   = 1.0,
            fwhm  = psf_fwhm,
            ratio = psf_ratio,
            x=x+eps(x),y=y+eps(y)
        )
        refs_M[:,i] .= psf_model.(
            rs[region_M] .* cos.(θs[region_M]),
            rs[region_M] .* sin.(θs[region_M]),
            # Eps to work around this bug: https://github.com/JuliaAstro/PSFModels.jl/issues/14
        ) .* (1 .- refs_are_rdi[i])
    end

    # f,ax,pl = Makie.scatterlines(
    #     phot_of_ref_at_target_inject_location,
    # )
    # Makie.hlines!(ax,[0.1,phot_of_target_at_inject_location])
    # display(f)

    # Include all frames with flux contamination in SNR-OPT, but feed frames with 
    # low contamination through a truncated SVD.
    frames_to_svd = phot_of_ref_at_target_inject_location .< 0.15phot_of_target_at_inject_location

    # Indices of "reference" images that have a model photometry equal
    # to the "target". Note that the target is actually included in the references
    # and multiple targets can be there at once. 
    # The "target" is really just used to set up some of the PSF modelling parameters.
    refs_with_maxphot = isapprox.(phot_of_ref_at_target_inject_location, phot_of_target_at_inject_location, atol=1e-4)
    println("# targets: ", count(refs_with_maxphot))
    if count(refs_with_maxphot) == 0
        @error "No target frames identified: something went wrong with the model"
        @show phot_of_target_at_inject_location
        @show phot_of_ref_at_target_inject_location
        return
    end

    # Do a intial subtraction, that is just the average of the "target" frames.
    # This will be our baseline.
    c = zeros(Float32, length(phot_of_ref_at_target_inject_location))
    c[refs_with_maxphot] .= 1.0
    phot_out = dot(c,phot_of_ref_at_target_inject_location)
    throughput = phot_out/phot_of_target_at_inject_location
    if !isfinite(throughput)
        @error "non finite throughput"
        return
    end
    if throughput ≈ 0
        @error "zero throughput"
        return
    end
    if any(!isfinite, c)
        @error "non finite coefficients"
        return
    end
    processed_S = refs_S*c
    processed_M = refs_M*c
    processed_S ./= throughput
    processed_M ./= throughput

    # Set the output to this initial averaged value. 
    # If we never improve, it will won't get overritten.
    # This allows us to do no subtraction at all if it doesn't help.
    out[region_S] .= processed_S
    fm[region_M] .= processed_M
    # display(imview(target))
    # display(imview(out))
    # display(imview(fm))
    # display(imview(fm .+ region_S))

    setglobal!(Main, :refcube,refcube)
    if length(filter(isfinite, processed_S)) == 0
        setglobal!(Main, :processed_S,processed_S)
        @warn "no pixels in subtraction region"
        return
    end
    t = collect(StatsBase.trim(filter(isfinite, processed_S);prop=0.05))
    if length(t) == 0
        @warn "no pixels in subtraction region"
        return
    end
    clip_std = std(t);  # Find std, while ignoring 10% most deviated pixels
    region_S_valid = all(
        (-4clip_std .< processed_S .< 4clip_std),
        dims=2
    )[:]
    # Don't do this clipping if we end up with nothing left, or take too much
    if count(region_S_valid) == 0 || mean(region_S_valid) < 0.25
        region_S_valid .= true
    end

    init_SNR = best_SNR = 1.0 ./ sqrt(mean(processed_S[region_S_valid].^2))

    # Consider whether we can bin frames before SNR-Opt, etc.
    # The SVD will handle this for frames with low contamination, but we can't apply
    # an SVD inteligently to create non-linear "SNR eigen modes". Instead we will just
    # bin down the data a bit.
    # This is especially usefull for multi-target SNR optimization where we might 
    # have many very similar reference images.
    # If the frames have < 2px of planet movement, and high correlation in the
    # optimization region, we can bin them together.
    # These will be a vector of vector of frames, grouped together by correlation and
    # angle.
    refs_to_bin = Tuple{Float64,Vector{Int}}[]
    for i in findall(.!frames_to_svd)
        # Bin data that is within 2% the same model flux and is well correlated
        model_flux_round = round(phot_of_ref_at_target_inject_location[i]*50)
        # If either there is nothing at this angle, or the first frame at that
        # angle is not very correlated, start a new group.
        found_matching_group = false
        for (model_flux′, group) in refs_to_bin
            if model_flux_round != model_flux′
                continue
            end
            # Compare correlation of current frame to first frame in the group
            j = group[1]
            c = cor(refs_O[:,i], refs_O[:,j])
            if c >= 0.99
                # println("highly correlated $i $j $c")
                push!(group, i)
                found_matching_group = true
                break
            end
        end
        if !found_matching_group
            push!(refs_to_bin, (model_flux_round, [i]))
        end
    end
    println("$(count(.!frames_to_svd)) frames with flux contamination binned to $(length(refs_to_bin))")

    # Now actually bin them
    refs_O_bin = zeros(eltype(refs_O), size(refs_O,1), length(refs_to_bin))
    refs_S_bin = zeros(eltype(refs_S), size(refs_S,1), length(refs_to_bin))
    phot_S_bin = zeros(eltype(phot_of_ref_at_target_inject_location), length(refs_to_bin))
    refs_M_bin = zeros(eltype(refs_M), size(refs_M,1), length(refs_to_bin))
    for (i,(_,group)) in enumerate(refs_to_bin)
        len = 0
        # println(phot_of_ref_at_target_inject_location[group])
        for j in group
            len += 1
            refs_O_bin[:,i] .+= @view refs_O[:,j]
            refs_S_bin[:,i] .+= @view refs_S[:,j]
            refs_M_bin[:,i] .+= @view refs_M[:,j]
            phot_S_bin[i] += phot_of_ref_at_target_inject_location[j]
        end
        # This was to average the different frames instead of summing them, but really there is no need to.
        # refs_O_bin[:,i] ./= len
        # refs_S_bin[:,i] ./= len
        # refs_M_bin[i] /= len
        # phot_S_bin[i] /= len
    end 

    # The truncated SVD takes a moderate amount of time.
    # We therefore do it only once, for the largest value of N_SVD and 
    # then index into it to select the right number of modes.
    # Note: if N_SVD=0 then we include everything. No need to SVD at all.
    if maximum(N_SVD) > 0 && count(frames_to_svd) > 0
        U, s, V = tsvd(Float32.(refs_O[:,frames_to_svd]), maximum(N_SVD))
    end

    # Loop through multiple N_SVD values and pick the best one.
    println("runtime         parameters      : SNR")
    println("\t\tinitial\t\t: $(best_SNR)")
    # Track the total time this region takes to process.
    t_tot = 0
    for n_SVD in N_SVD
        # Track the time of this block
        t = @elapsed begin
            # If we are using all frames, we do not use the SVD.
            if n_SVD == 0 || count(frames_to_svd) == 0
                refs_O_prepared = refs_O_bin
                refs_S_prepared = refs_S_bin 
                refs_M_prepared = refs_M_bin
                phot_S_prepared = phot_S_bin
            else
                # These are the N largest SVD components among the frames that do not have significant
                # flux overlap with the planet PSF (<10%).
                # U gives the reference region eigen-images.
                # We then "rotate" (or "transform") the modelled photometry and subtraction
                # regions to match the new truncated eigen-image basis.
                # We have to multiply by the inverse of the singular value vector to maintain
                # the correct intensity scales between these three products.
                refs_O_low_contam_rotated = U[:,1:n_SVD]
                refs_S_low_contam_rotated = refs_S[:,frames_to_svd] * V[:,1:n_SVD] * inv(Diagonal(s[1:n_SVD]))
                refs_M_low_contam_rotated = refs_M[:,frames_to_svd] * V[:,1:n_SVD] * inv(Diagonal(s[1:n_SVD]))
                phot_S_low_contam_rotated = (phot_of_ref_at_target_inject_location[frames_to_svd]'*( V[:,1:n_SVD] * inv(Diagonal(s[1:n_SVD]))))'

                # Append SVD'd frames to regular contaminated frames.
                # We will then proceed with SNR-Opt using this hybrid basis of SVD eigen-images and
                # natural frames that contain significant flux overlap.
                refs_O_prepared = Float32[refs_O_low_contam_rotated;; refs_O_bin] #refs_O[:,.!frames_to_svd]]
                refs_S_prepared = Float32[refs_S_low_contam_rotated;; refs_S_bin] #refs_S[:,.!frames_to_svd]]
                refs_M_prepared = Float32[refs_M_low_contam_rotated;; refs_M_bin] #refs_S[:,.!frames_to_svd]]
                phot_S_prepared = Float32[phot_S_low_contam_rotated;  phot_S_bin] #phot_of_ref_at_target_inject_location[.!frames_to_svd]]
            end

            # refs_O_rotated = refs_O[:,:]
            # refs_S_rotated = refs_S
            # phot_S_rotated = phot_of_ref_at_target_inject_location

            # Maximizing or minimizing the SNR is actually equivalent, as long 
            # as we do throughput correction! One makes the planet positive,
            # the other makes the planet negative.
            (f, f′) = make_snr_callbacks(refs_O_prepared, phot_S_prepared)
            # (f, f′) = make_nsr_callbacks(refs_O_prepared, phot_S_prepared)

            # Start with random coefficients, but seed with the shape of the region
            # for some reprodicibility.
            Random.seed!(hash(region_S))
            c = rand(size(refs_S_prepared,2))

            # In-place gradient function for optimizer
            g!(G, c) = G .= f′(c)
            # Run optimization
            try
                result = Optim.optimize(
                    f,
                    g!,
                    convert(Vector{Float32}, c),
                    LBFGS(),
                    Optim.Options(
                        show_trace=false,
                        show_every=1000,
                        iterations = 30000,
                        g_reltol = 1e-12,                
                        g_abstol = 1e-12,        
                        f_reltol = 1e-12,                
                        f_abstol = 1e-12,                
                    ),
                )

                # Even if the result is not converged it is often pretty good. 
                # Emit a warning, but still go on to check if these coefficients
                # have improved the reduction.
                if !Optim.converged(result)
                    @warn "region did not converge"
                    display(result)
                end
                c = Optim.minimizer(result)
            catch err
                @error "Error during optimization" exception=(err, catch_backtrace())
            end

            # Normal output
            processed_S = refs_S_prepared*c
            # Model output:
            processed_M = refs_M_prepared*c

            phot_out = dot(c, phot_S_prepared)
            
            throughput = phot_out/phot_of_target_at_inject_location

            # Throughput normalize (just once for full anulus)
            # TODO: calculate throughput at inner edge, outer edge, and middle,
            # then do linear interpolation. That reduced the edge artifacts.
            processed_S ./= throughput
            processed_M ./= throughput

            clip_std = std(StatsBase.trim(vec(processed_S); prop=0.05)) # Find std, while ignoring 10% most deviated pixels
            region_S_valid = all(
                (-4clip_std .< processed_S .< 4clip_std),
                dims=2
            )[:]
            # Don't do this clipping if we end up with nothing left, or take too much
            if count(region_S_valid) == 0 || mean(region_S_valid) < 0.25
                region_S_valid .= true
            end

            new_SNR = 1.0 ./ sqrt(mean(processed_S[region_S_valid].^2))
        end
        t_tot += t
        print("$t\tn_SVD=$n_SVD \t: $(new_SNR)")
        if new_SNR > best_SNR
            print("\t🟢")
            best_SNR = new_SNR
            out[region_S] .= processed_S
            fm[region_M] .= processed_M
        else
            print("\t🔻")
        end
        println()

        if n_SVD >= count(frames_to_svd)
            break
        end
    end
    improvement = best_SNR/init_SNR
    print("$t_tot\timprovement \t: ×$improvement")
    if improvement > 100
        println("\t🤩")
    elseif improvement > 50
        println("\t😄")
    elseif improvement > 10
        println("\t🙂")
    elseif improvement > 1
        println("\t🤷")
    elseif improvement < 1.5
        println("\t😒")
    end
    println()
end


function _snr(
    c::AbstractVector,
    imgdat,
    photdat,
    temp2,
    temp3,
)
    phot = dot(c, photdat)
    mul!(temp2, imgdat, c)
    temp3 .= temp2.^2
    noise = sqrt(mean(temp3))
    # There is a free parameter for SNR, where all values can be scaled up or down by a constant 
    # factor. This can make the optimizer run away to large eg 1e9 values.
    # Add a very slight L2 regularization penalty on the coefficients to
    # keep their values closer to zero and prevent numerical issues.
    return phot/noise + 0.001dot(c,c)
end

function _nsr(
    c::AbstractVector,
    imgdat,
    photdat,
    temp2,
    temp3,
)
    phot = dot(c, photdat)
    mul!(temp2, imgdat, c)
    temp3 .= temp2.^2
    noise = sqrt(mean(temp3))
    return noise/phot
end
"""
    (snr, snr′,) = make_snr_callbacks(refs_O::Matrix, phot_S::Vector)

Given image data and model photometry, return a callback for calculating
the SNR given linear combination coefficients and a callback for calculating
the gradient.

"""
function make_snr_callbacks(refs_O_rotated::AbstractMatrix{<:Number}, phot_S_rotated::AbstractVector{<:Number})
    imgdat  = convert(Matrix{Float32}, refs_O_rotated)
    photdat = convert(Vector{Float32}, phot_S_rotated)
    temp2=similar(imgdat, size(imgdat,1))
    temp3=similar(imgdat, size(imgdat,1))

    function snr(
        c::AbstractVector,
        imgdat=imgdat,
        photdat=photdat,
        temp2=temp2,
        temp3=temp3,
    )
        return _snr(c, imgdat, photdat, temp2, temp3)
    end
    
    # In previous versions of the SNAP pipeline, we used a convex optimization library
    # to formulate this problem. In practice, a simple iterative gradient optimization 
    # seems to perform pretty well too.

    # We use the Enzyme library to automatically calculate the gradient of the SNR 
    # with respect to the coefficients provided.
    # A hand-coded gradient calculation might be faster. I have the math worked out
    # somewhere for the analytic gradient, but would have to dig it out from the 
    # unknown depths of an old laptop hard drive, or just re-calculate it by hand.

    # Allocate temporary storage used by the SNR calculation and Enzyme gradient
    # calculation once, and re-use for all calls into this function.
    # Note that it is not thread-safe as currently implemented!! (would need to use
    # task-local storage for these scratch spaces.)
    temp2_diff = copy(temp2)
    temp3_diff = copy(temp3)
    temp2_d = Duplicated(temp2, temp2_diff)
    temp3_d = Duplicated(temp3, temp3_diff)
    c_diff = zeros(Float32, length(phot_S_rotated))
    function snr′(c::AbstractVector{Float32})
        fill!(temp2, 0)
        fill!(temp3, 0)
        fill!(temp2_diff, 0)
        fill!(temp3_diff, 0)
        fill!(c_diff, 0)
        # Enzyme can also calculate both the primal value and the gradient 
        # simultaneously, however, last I checked it was not very fast.
        out = Enzyme.autodiff(
            Reverse,
            _snr,
            Duplicated(c, c_diff), 
            Const(imgdat),
            Const(photdat),
            temp2_d,
            temp3_d,
        )
        return c_diff
    end
    (snr, snr′)
end



function make_nsr_callbacks(refs_O_rotated::AbstractMatrix{<:Number}, phot_S_rotated::AbstractVector{<:Number})
    imgdat  = convert(Matrix{Float32}, refs_O_rotated)
    photdat = convert(Vector{Float32}, phot_S_rotated)
    temp2=similar(imgdat, size(imgdat,1))
    temp3=similar(imgdat, size(imgdat,1))

    function nsr(
        c::AbstractVector,
        imgdat=imgdat,
        photdat=photdat,
        temp2=temp2,
        temp3=temp3,
    )
        return _nsr(c, imgdat, photdat, temp2, temp3)
    end
    temp2_diff = copy(temp2)
    temp3_diff = copy(temp3)
    temp2_d = Duplicated(temp2, temp2_diff)
    temp3_d = Duplicated(temp3, temp3_diff)
    c_diff = zeros(Float32, length(phot_S_rotated))
    function nsr′(c::AbstractVector{Float32})
        fill!(temp2, 0)
        fill!(temp3, 0)
        fill!(temp2_diff, 0)
        fill!(temp3_diff, 0)
        fill!(c_diff, 0)
        out = Enzyme.autodiff(
            Reverse,
            _nsr,
            Duplicated(c, c_diff), 
            Const(imgdat),
            Const(photdat),
            temp2_d,
            temp3_d,
        )
        return c_diff
    end
    (nsr, nsr′)
end