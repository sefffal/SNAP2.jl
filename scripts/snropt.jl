using Random: Random
using Enzyme
using Optim: Optim


function snropt_all(fnames_pattern, rotthreshpx, regions_S, regions_O; force=false, kwargs...)
    fnames = glob(fnames_pattern)
    imgs = load.(fnames)
    refcube=stack(imgs)
    i = 0
    for (fname, targ) in zip(fnames, imgs)
        i+=1
        outfname = replace(fname, ".fits"=>".snropt.fits")
        if !force && isfile(outfname) &&  Base.Filesystem.mtime(outfname) > Base.Filesystem.mtime(fname)
            out = load(outfname)
        else
            out = deepcopy(targ)
            fill!(out, NaN)
        end

        # for (reg_S, reg_O) in zip(regions_S, regions_O)
        # Threads.@threads :dynamic 
        for (reg_S, reg_O) in collect(zip(regions_S, regions_O))
            snropt_region!(out, targ, refs, reg_S.>0, reg_O.>0, refcube; kwargs...)
        end
        AstroImages.writefits(outfname, out)      
        println(outfname, "\t($i)")  
    end
end



function snropt_frame(fname, refnames_pattern, region_S, region_O;force=false, kwargs...)
    target = load(fname)
    outfname = replace(fname, ".fits"=>".snropt.fits")
    if !force && isfile(outfname) &&  Base.Filesystem.mtime(outfname) > Base.Filesystem.mtime(fname)
        out = load(outfname)
    else
        out = deepcopy(target)
        fill!(out, NaN)
    end
    refnames = glob(refnames_pattern)
    refs = load.(refnames)
    refcube=stack(refs)
    i = 0
    for (reg_S, reg_O) in zip(region_S, region_O)
        i+= 1
        println("region #$i")
        snropt_region!(out, target, refs, reg_S.>0, reg_O.>0, refcube; kwargs...)
    end
    AstroImages.writefits(outfname, out)
    return out
end



function snropt_region!(
        out::AstroImage,
        target::AstroImage,
        refs::AbstractVector{<:AstroImage},
        region_S::AbstractMatrix{Bool},
        region_O::AbstractMatrix{Bool},
        refcube=stack(refs);
        N_SVD=0,
        psf_fwhm,
        psf_ratio=0.0,
)
    rs = imgsep(target)
    θs = imgang(target)
    # sep = mean(rs[region_S])
    # angle = deg2rad(target["ANGLE_MEAN"]::Float64)
    # angles = deg2rad.(getindex.(refs, "ANGLE_MEAN")::Vector{Float64})
    # valid_II = findall(allowed)

    ## Given a list of frames, calculate the photometry of the planet at the mean location
    # of the target frame.
    # We can use that for throughput normalization, for tuning hyper-parameters, and for 
    # SNAP (down the road)
    # Issue: we really need to be able to rotate this into the SVD basis.

    # We simulate the planet at the average PA and average separation of the target subtraction
    # region.
    inject_planet_sep_ave = mean(rs[region_S])
    inject_planet_pa_ave = deg2rad.(Snap.meandegrees(vec(rad2deg.(θs[region_S]))))
    # We now figure out how much the planet has rotated in each reference frame.
    # TODO: confirm the sign is correct in second part of expression (no effect if PSF is symmetric through.)
    ref_planet_pa = inject_planet_pa_ave .+ deg2rad.(getindex.(refs, "ANGLE_MEAN").- target["ANGLE_MEAN"] )
    ref_planet_pa .= rem2pi.(ref_planet_pa, RoundDown)

    # for now don't consider scaling due to SDI, only consider rotation from ADI.
    x = inject_planet_sep_ave * cos(inject_planet_pa_ave)
    y = inject_planet_sep_ave * sin(inject_planet_pa_ave)
    psf_model = PSFModels.airydisk(;
        amp   = 1.0,
        fwhm  = psf_fwhm,
        ratio = psf_ratio,
        x,y
    )

    # Calculate just the point photometry at the peak of the target
    phot_of_target_at_inject_location = psf_model(x,y)
    phot_of_ref_at_target_inject_location = psf_model.(
        inject_planet_sep_ave .* cos.(ref_planet_pa),
        inject_planet_sep_ave .* sin.(ref_planet_pa),
    )

    # f,ax,pl = Makie.scatterlines(
    #     ref_planet_pa
    # )
    # Makie.hlines!(ax,inject_planet_pa_ave)
    # display(f)

    # Calculate actual models to see what the PSF looks like
    # targ_M = psf_model.(
    #     rs[region_S] .* cos.(θs[region_S]),
    #     rs[region_S] .* sin.(θs[region_S]),
    # )
    # refs_M = zeros(eltype(refcube), count(region_S), length(valid_II))
    # for i in eachindex(valid_II)
    #     x = inject_planet_sep_ave * cos(ref_planet_pa[i])
    #     y = inject_planet_sep_ave * sin(ref_planet_pa[i])
    #     psf_model = PSFModels.airydisk(;
    #         amp   = 1.0,
    #         fwhm  = psf_fwhm,
    #         ratio = psf_ratio,
    #         x,y
    #     )
    #     refs_M[:,i] .= psf_model.(
    #         rs[region_S] .* cos.(θs[region_S]),
    #         rs[region_S] .* sin.(θs[region_S]),
    #     )
    # end
    # outs = map(axes(refs_M,2)) do i
    #     o = deepcopy(out)
    #     o[region_S] .= refs_M[:,i]
    #     return o
    # end
    # out[region_S] .= targ_M
    # return out, outs


    targ_O = target[region_O]
    refs_O = zeros(eltype(refcube), count(region_O), length(refs))
    refs_S = zeros(eltype(refcube), count(region_S), length(refs))
    for i in eachindex(refs)
        refs_O[:,i] .= refcube[region_O,i]
        refs_S[:,i] .= refcube[region_S,i]
    end
    valid = vec(all(isfinite,refs_O,dims=2)) .& isfinite.(targ_O)

    # Reject pixels that are strongly deviated from the clipped standard deviation.
    clip_std = std(StatsBase.trim(vec(view(refs_O,valid,:)); prop=0.05)) # Find std, while ignoring 5% most deviated pixels
    valid .&= all(
        (-5clip_std .< view(refs_O,valid,:) .< 5clip_std) .&
        (-5clip_std .< view(targ_O,valid) .< 5clip_std),
        dims=2
    ) 
    # Don't do this clipping if we end up with nothing left, or take too much
    if count(valid) == 0|| mean(valid) < 0.25
        valid .= true
    end




    # Identify the pixels that are safe to use for evaluating the SNR of the subtraction
    # Avoid looking at bad pixels
    # Ignore the bottom and top 0.5% of pixels
    clip_ex = extrema(StatsBase.trim(vec(refs_S); prop=0.01)) # Find std, while ignoring 5% most deviated pixels
    region_S_valid = all(
        (clip_ex[1] .< refs_S .< clip_ex[2]) .&
        (clip_ex[1] .< target[region_S] .< clip_ex[2]),
        dims=2
    )[:]
    # Don't do this clipping if we end up with nothing left, or take too much
    if count(region_S_valid) == 0 || mean(region_S_valid) < 0.25
        region_S_valid .= true
    end

    f,ax,pl = Makie.scatterlines(
        phot_of_ref_at_target_inject_location,
    )
    Makie.hlines!(ax,[0.1,phot_of_target_at_inject_location])
    display(f)


    # Include all frames with flux contamination in SNR-OPT, but feed frames with 
    # low contamination through a truncated SVD.
    frames_to_svd = phot_of_ref_at_target_inject_location .< 0.1
    
    # Loop through multiple N_SVD values and pick the best one.
    # We start with no subtraction at all.
    out[region_S] .= target[region_S]
    init_SNR = best_SNR = 1.0 / sqrt(mean(target[region_S][region_S_valid].^2))
    println("initial\t\t: $(best_SNR)\t*")
    for n_SVD in N_SVD
        if n_SVD == 0 || count(frames_to_svd) == 0
            refs_O_rotated = refs_O[valid,:]
            refs_S_rotated = refs_S
            phot_S_rotated = phot_of_ref_at_target_inject_location
        else
            U, s, V = tsvd(Float32.(refs_O[valid,frames_to_svd]), n_SVD)
            refs_O_low_contam_rotated = U
            refs_S_low_contam_rotated = refs_S[:,frames_to_svd] *V * inv(Diagonal(s))
            phot_S_low_contam_rotated = (phot_of_ref_at_target_inject_location[frames_to_svd]'*(V*inv(Diagonal(s))))'

            # Append SVD'd frames to regular contaminated frames
            refs_O_rotated = Float32[refs_O_low_contam_rotated;; refs_O[valid,.!frames_to_svd]]
            refs_S_rotated = Float32[refs_S_low_contam_rotated;; refs_S[:,.!frames_to_svd]]
            phot_S_rotated = Float32[phot_S_low_contam_rotated; phot_of_ref_at_target_inject_location[.!frames_to_svd]]
        end

        # refs_O_rotated = refs_O[valid,:]
        # refs_S_rotated = refs_S
        # phot_S_rotated = phot_of_ref_at_target_inject_location

        # (snr, snr′) = make_snr_callbacks(refs_O_rotated, phot_S_rotated)
        (f, f′) = make_nsr_callbacks(refs_O_rotated, phot_S_rotated)

        # return (refs_O_rotated, phot_S_rotated, targ_O[valid], f, f′)

        # Start with random coefficients, but seed with the shape of the region
        # for some reprodicibility
        # rng = Random.Xoshiro(hash(region_S))
        rng = Random.Xoshiro(count(region_S))
        c = randn(rng, size(refs_S_rotated,2))
        g!(G, c) = G .= f′(c)
        result = Optim.optimize(
            f,
            g!,
            convert(Vector{Float32}, c),
            LBFGS(),
            Optim.Options(show_trace=false, show_every=1000, iterations = 15000)
        )

        if !Optim.converged(result)
            @warn "region did not converge"
            display(result)
        end
        c = Optim.minimizer(result)


        # Normal output
        processed_S = target[region_S] .- refs_S_rotated*c

        # Model output:
        # processed_S .= targ_M .- vec(sum(reshape(c,1,:) .* refs_M_rotated,dims=2))

        # display(imview(out))
        phot_out = dot(c,phot_S_rotated)
        
        throughput = phot_out/phot_of_target_at_inject_location

        # Throughput normalize (just once for full anulus)
        processed_S ./= throughput

        new_SNR = 1.0 ./ sqrt(mean(processed_S[region_S_valid].^2))
        print("n_SVD=$n_SVD \t: $(new_SNR)")
        if new_SNR > best_SNR
            print("\t*")
            best_SNR = new_SNR
            out[region_S] .= processed_S
        end
        println()
    end
    println("improvement = $(best_SNR/init_SNR)")
    println()

    # @show throughput phot_out phot_of_target_at_inject_location
    # display(
    #     Makie.scatterlines(phot_of_ref_at_target_inject_location)
    # )

    # return out .* region_S
    # return c

    # @show size(valid) count(valid) size(refs_O_rotated)


    # modes = zeros(size(region_O)...,size(refs_O_rotated,2))
    # modes_region_O = @view modes[region_O,:]
    # modes_region_O[valid,:] .=refs_O_rotated
    # modes_region_S = @view modes[region_S,:]
    # modes_region_S .= refs_S_rotated
    
    # return modes

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
    return phot/noise
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
    return phot/noise
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


    # function snr′_with_primal(c::AbstractVector{Float32})
    #     fill!(temp1, 0)
    #     fill!(temp2, 0)
    #     fill!(temp3, 0)
    #     fill!(temp1_diff, 0)
    #     fill!(temp2_diff, 0)
    #     fill!(temp3_diff, 0)
    #     fill!(c_diff, 0)
    #     out, primal = Enzyme.autodiff(
    #         ReverseWithPrimal,
    #         snr,
    #         Duplicated(c, c_diff), 
    #         Const(imgdat),
    #         Const(photdat),
    #         temp1_d,
    #         temp2_d,
    #         temp3_d,
    #     )
    #     return primal::eltype(c), c_diff
    # end


    # optimize(f, g!, x0, LBFGS())

    # oh = Enzyme.onehot(c_diff)
    # function snr′(c::AbstractVector{Float32})
    #     fill!(temp1, 0)
    #     fill!(temp2, 0)
    #     fill!(temp3, 0)
    #     fill!(temp1_diff, 0)
    #     fill!(temp2_diff, 0)
    #     fill!(temp3_diff, 0)
    #     primal, out = Enzyme.autodiff(
    #         Enzyme.Forward,
    #         snr,
    #         Enzyme.BatchDuplicated,
    #         Enzyme.BatchDuplicated(c,oh),
    #         Const(imgdat),
    #         Const(photdat),
    #         temp1_dnn,
    #         temp2_dnn,
    #         temp3_dnn
    #     )
    #     c_diff .= tuple(out...)
    #     return primal, c_diff
    # end

        
    (snr, snr′)#, snr′_with_primal)
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


    # c_diff = zeros(eltype(photdat), length(photdat))
    # temp2_diff = copy(temp2)
    # temp3_diff = copy(temp3)
    # temp2_d = Duplicated(temp2, temp2_diff)
    # temp3_d = Duplicated(temp3, temp3_diff)
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


    (nsr, nsr′)#, snr′_with_primal)
end