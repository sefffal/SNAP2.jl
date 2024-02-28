using FITSIO
using StatsBase
using TSVD # truncated SVD library
using PSFModels: PSFModels
function loci2_region(fname, refnames_pattern::AbstractString, rotthreshpx, region_S, region_O; kwargs...)
    target = load(fname)
    outfname = replace(fname, ".fits"=>".sub.fits")
    if isfile(outfname)
        out = load(outfname)
    else
        out = deepcopy(target)
        fill!(out, NaN)
    end
    refnames = glob(refnames_pattern)
    refs = load.(refnames)
    out = loci2_region!(out, target, refs, rotthreshpx, region_S, region_O; kwargs...)
    AstroImages.writefits(outfname, out)
end
using LinearSolve
function loci2_region!(
    out::AstroImage,
    target::AstroImage,
    refs::AbstractVector{<:AstroImage},
    rotthreshpx::Number,
    region_S::AbstractMatrix{Bool},
    region_O::AbstractMatrix{Bool},
    refcube=stack(refs);
    N_SVD=60,
    psf_fwhm,
    psf_ratio=0.0,
)
    rs = imgsep(target)
    θs = imgang(target)
    sep = mean(rs[region_S])
    angle = deg2rad(target["ANGLE_MEAN"]::Float64)
    angles = deg2rad.(getindex.(refs, "ANGLE_MEAN")::Vector{Float64})
    allowed = sep .* rem2pi.(abs.(angle .- angles), RoundDown) .> rotthreshpx
    valid_II = findall(allowed)

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
    ref_planet_pa = inject_planet_pa_ave .+ deg2rad.(getindex.(refs[valid_II], "ANGLE_MEAN").- target["ANGLE_MEAN"] )
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
    refs_O = zeros(eltype(refcube), count(region_O), length(valid_II))
    refs_S = zeros(eltype(refcube), count(region_S), length(valid_II))
    for (i,j) in enumerate(valid_II)
        refs_O[:,i] .= refcube[region_O,j]
        refs_S[:,i] .= refcube[region_S,j]
    end
    valid = vec(all(isfinite,refs_O,dims=2)) .& isfinite.(targ_O)

    # Reject pixels that are strongly deviated from the clipped standard deviation.
    clip_std = std(StatsBase.trim(vec(view(refs_O,valid,:)); prop=0.05)) # Find std, while ignoring 5% most deviated pixels
    valid .&= all((-5clip_std .< view(refs_O,valid,:) .< 5clip_std) .& (-5clip_std .< view(targ_O,valid) .< 5clip_std),dims=2) 

    if count(valid) == 0
        out[region_S] .= target[region_S]
        return out
    end

    ## SVD
    # U, s, V = tsvd(Float32.(refs_O[valid,:]), N_SVD)
    # # U, s, V = svd(Float32.(refs_O[valid,:]))
    # refs_O_rotated = U
    # @show size(U) size(V) size(s) size(refs_S)
    # refs_S_rotated = refs_S *V * inv(Diagonal(s))
    # phot_S_rotated = phot_of_ref_at_target_inject_location[valid_II]'*(V*inv(Diagonal(s)))
    # # phot_S_rotated = (V'*inv(Diagonal(s)))*phot_of_ref_at_target_inject_location[valid_II]'
    # refs_M_rotated =  refs_M*V*inv(Diagonal(s))


    refs_O_rotated = refs_O[valid,:]
    refs_S_rotated = refs_S
    # refs_M_rotated = refs_M
    phot_S_rotated = phot_of_ref_at_target_inject_location




    # U, s, V = tsvd(Float32.([refs_O[valid,:];  refs_S]), N_SVD)
    # refs_O_rotated = U[1:count(valid),:]
    # refs_S_rotated = U[count(valid)+1:end,:]
    # phot_S_rotated = phot_of_ref_at_target_inject_location[valid_II]'*V

    # U is now our reference optimization zone.
    # We also want to rotate our subtraction zone in the same way, so that we will be able to apply the correction.
    # refs_S_rotated = refs_S*V

    ## Linear solve
    c = refs_O_rotated \ targ_O[valid]
    # prob = LinearProblem(refs_O_rotated, Float32.(targ_O[valid]))
    # @showtime sol = solve(prob, SVDFactorization())#, KrylovJL_GMRES())
    # sol = solve(prob, QRFactorization())#, KrylovJL_GMRES())
    # @showtime sol = solve(prob, NormalCholeskyFactorization())#, KrylovJL_GMRES())
    # @showtime sol = solve(prob, RFLUFactorization())#, KrylovJL_GMRES())
    # @showtime sol = solve(prob, AppleAccelerateLUFactorization())#, KrylovJL_GMRES())
    # @showtime sol = solve(prob, SimpleLUFactorization())#, KrylovJL_GMRES())
    # c = sol.u

    # ## Linear solve
    # # c = refs_O[valid,:] \ targ_O[valid]
    # prob = LinearProblem(Float32.(refs_O[valid,:]), Float32.(targ_O[valid]))
    # # @showtime sol = solve(prob, SVDFactorization())#, KrylovJL_GMRES())
    # sol = solve(prob, QRFactorization())#, KrylovJL_GMRES())
    # # @showtime sol = solve(prob, NormalCholeskyFactorization())#, KrylovJL_GMRES())
    # # @showtime sol = solve(prob, RFLUFactorization())#, KrylovJL_GMRES())
    # # @showtime sol = solve(prob, AppleAccelerateLUFactorization())#, KrylovJL_GMRES())
    # # @showtime sol = solve(prob, SimpleLUFactorization())#, KrylovJL_GMRES())
    # c = sol.u



    # Normal output
    out[region_S] .= target[region_S] .- vec(sum(reshape(c,1,:) .* refs_S_rotated,dims=2))

    # Model output:
    # out[region_S] .= targ_M .- vec(sum(reshape(c,1,:) .* refs_M_rotated,dims=2))

    # display(imview(out))
    phot_out = phot_of_target_at_inject_location .- sum(c .* phot_S_rotated)
    
    throughput = phot_out/phot_of_target_at_inject_location

    # Throughput normalize (just once for full anulus)
    out[region_S] ./= throughput

    # @show throughput phot_out phot_of_target_at_inject_location
    # display(
    #     Makie.scatterlines(phot_of_ref_at_target_inject_location)
    # )

    return out .* region_S

    # @show size(valid) count(valid) size(refs_O_rotated)


    # modes = zeros(size(region_O)...,size(refs_O_rotated,2))
    # modes_region_O = @view modes[region_O,:]
    # modes_region_O[valid,:] .=refs_O_rotated
    # modes_region_S = @view modes[region_S,:]
    # modes_region_S .= refs_S_rotated
    
    # return modes

end

function loci2_frame(fname, refnames_pattern, rotthreshpx, region_S, region_O; kwargs...)
    target = load(fname)
    outfname = replace(fname, ".fits"=>".sub.fits")
    if isfile(outfname)
        out = load(outfname)
    else
        out = deepcopy(target)
        fill!(out, NaN)
    end
    refnames = glob(refnames_pattern)
    refs = load.(refnames)
    refcube=stack(refs)
    for (reg_S, reg_O) in zip(region_S, region_O)
        loci2_region!(out, target, refs, rotthreshpx, reg_S.>0, reg_O.>0, refcube; kwargs...)
    end
    AstroImages.writefits(outfname, out)
    return out
end



function loci2_all(fnames_pattern, rotthreshpx, regions_S, regions_O; force=false, kwargs...)
    fnames = glob(fnames_pattern)
    imgs = load.(fnames)
    refcube=stack(imgs)
    i = 0
    for (fname, targ) in zip(fnames, imgs)
        i+=1
        ii = setdiff(axes(refcube,3),i)
        refs_this = @view imgs[ii]
        refcube_this = @view refcube[:,:,ii]

        outfname = replace(fname, ".fits"=>".sub.fits")
        if !force && isfile(outfname) &&  Base.Filesystem.mtime(outfname) > Base.Filesystem.mtime(fname)
            out = load(outfname)
        else
            out = deepcopy(targ)
            fill!(out, NaN)
        end

        # for (reg_S, reg_O) in zip(regions_S, regions_O)
        # Threads.@threads :dynamic 
        for (reg_S, reg_O) in collect(zip(regions_S, regions_O))
            loci2_region!(out, imgs[i], refs_this, rotthreshpx, reg_S.>0, reg_O.>0, refcube_this; kwargs...)
        end
        AstroImages.writefits(outfname, out)      
        println(outfname, "\t($i)")  
    end
end