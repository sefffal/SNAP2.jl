using FITSIO
using StatsBase
using TSVD # truncated SVD library
using PSFModels: PSFModels

export loci2_region, loci2_region!, loci2_frame, loci2_all

#=

PSF Modelling

We want to replace the synthetic PSF with the real unsaturated PSF template---or at least a more
realistic synthetic one that has asymetries.

Since LOCI is normally done in the speckle-aligned rotation frame, this is kind of easy here.
But in theory someone could apply LOCI to a dataset that has already been rotated.
This also is needed to generalize to multi-targ SNROpt where the speckles have rotated.

So we need to figure out if we should rotate the template. How can we do this?

ANGLE_MEAN gives the rotation we applied or need to apply to put the planets North up.
What about putting planets.

Taken to the max: what if we just create a PSF template for every file..?
We could apply all the same transformations to it; we could see what's happening.
We already have a model PSF per region coming *OUT* of SNRopt.
What if we put a model PSF per region going *IN* to subtractions?


=#

function loci2_region(fname, refnames_pattern::AbstractString, rotthreshpx, region_S, region_O; kwargs...)
    target = load(fname)
    outfname = replace(fname, ".fits"=>".loci.fits")
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
    rotthreshpx,
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
    angle = deg2rad(Float64(target["ANGLE_MEAN"])::Float64)
    angles = deg2rad.(Float64.(getindex.(refs, "ANGLE_MEAN"))::Vector{Float64})

    out[region_S] .= target[region_S]

    if count(region_S) == 0 || count(region_O) == 0
        @warn "cannot perform PSF subtraction with an empty region " count(region_S) count(region_O)
    end

    if length(filter(isfinite, target[region_S])) == 0
        @warn "empty or all-NaN subtraction region"
        return
    end

    # We simulate the planet at the average PA and average separation of the target subtraction
    # region.
    inject_planet_sep_ave = mean(rs[region_S])
    inject_planet_pa_ave = deg2rad.(SNAP.meandegrees(vec(rad2deg.(θs[region_S]))))
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

    targ_O = target[region_O]

   
    pixels = collect(StatsBase.trim(vec(target[region_S]); prop=0.05))
    if length(pixels) <= 0
        @warn "clipped all pixels"
        return
    end
    clip_std = std(pixels) # Find std, while ignoring 10% most deviated pixels
    region_S_valid = all(
        (-4clip_std .< target[region_S] .< 4clip_std),
        dims=2
    )[:]
    # Don't do this clipping if we end up with nothing left, or take too much
    if count(region_S_valid) == 0 || mean(region_S_valid) < 0.25
        region_S_valid .= true
    end

    # TODO: it would be better to work in fractions of psf_fwhm.
    # We also should de-duplicate. There are often frames where we are spinning
    # fast enough that a threshold of x and a threshold of y px work out to the 
    # same frames.

    # Loop through optimizing the rejection thrreshold ratio
    init_SNR = best_SNR = 1.0 / sqrt(mean(target[region_S][region_S_valid].^2))
    println("initial\t\t\t\t\t\t: $(best_SNR)\t*")
    for rotthreshpx_val in rotthreshpx
        distances =  sep .* tan.(rem2pi.(abs.(angle .- angles), RoundDown))
        allowed = distances .> rotthreshpx_val
        valid_II = findall(allowed)
        if isempty(valid_II)
            continue
        end

        ## Given a list of frames, calculate the photometry of the planet at the mean location
        # of the target frame.
        # We can use that for throughput normalization, for tuning hyper-parameters, and for 
        # SNAP (down the road)
        # Issue: we really need to be able to rotate this into the SVD basis.

        # We now figure out how much the planet has rotated in each reference frame.
        # TODO: confirm the sign is correct in second part of expression (no effect if PSF is symmetric through.)
        ref_planet_pa = inject_planet_pa_ave .+ deg2rad.(getindex.(refs[valid_II], "ANGLE_MEAN").- target["ANGLE_MEAN"] )
        ref_planet_pa .= rem2pi.(ref_planet_pa, RoundDown)
        phot_of_ref_at_target_inject_location = psf_model.(
            inject_planet_sep_ave .* cos.(ref_planet_pa),
            inject_planet_sep_ave .* sin.(ref_planet_pa),
        )

        # TODO: could hoist this out and then just do views into it
        refs_O = zeros(eltype(refcube), count(region_O), length(valid_II))
        refs_S = zeros(eltype(refcube), count(region_S), length(valid_II))
        for (i,j) in enumerate(valid_II)
            refs_O[:,i] .= refcube[region_O,j]
            refs_S[:,i] .= refcube[region_S,j]
        end
        valid = vec(all(isfinite,refs_O,dims=2)) .& isfinite.(targ_O)
        if count(valid) <= 0
            continue
        end
        # Reject pixels that are strongly deviated from the clipped standard deviation.
        clip_std = std(StatsBase.trim(vec(view(refs_O,valid,:)); prop=0.05)) # Find std, while ignoring 5% most deviated pixels
        valid .&= all((-5clip_std .< refs_O .< 5clip_std) .& (-5clip_std .< targ_O .< 5clip_std),dims=2) 
    
        if count(valid) <= 0
            continue
        end

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

        if any(!isfinite, refs_O[valid,:])
            error("invalid condition: non finite pixels in masked optimization region")
        end

        local U,s,V
        if maximum(N_SVD) > 0
            try
                U, s, V = tsvd(Float32.(refs_O[valid,:]), min(maximum(N_SVD), size(refs_O,2)))
            catch err
                @error "Issue during truncated SVD (TSVD.jl)" exception=(err, catch_backtrace())
                return
            end
        end

        for n_SVD in N_SVD
            if maximum(N_SVD) > 0 && n_SVD > size(U, 2)
                continue
            end

            ## SVD
            if n_SVD == 0
                refs_O_rotated = refs_O[valid,:]
                refs_S_rotated = refs_S 
                # refs_M_rotated = refs_M_bin
                phot_S_rotated = phot_of_ref_at_target_inject_location
            else
                refs_O_rotated = U[:,1:n_SVD]
                refs_S_rotated = refs_S[:,:] * V[:,1:n_SVD] * inv(Diagonal(s[1:n_SVD]))
                phot_S_rotated = (phot_of_ref_at_target_inject_location[:]'*( V[:,1:n_SVD] * inv(Diagonal(s[1:n_SVD]))))'
            end


            # refs_O_rotated = refs_O[valid,:]
            # refs_S_rotated = refs_S
            # # refs_M_rotated = refs_M
            # phot_S_rotated = phot_of_ref_at_target_inject_location




            # U, s, V = tsvd(Float32.([refs_O[valid,:];  refs_S]), n_SVD)
            # refs_O_rotated = U[1:count(valid),:]
            # refs_S_rotated = U[count(valid)+1:end,:]
            # phot_S_rotated = phot_of_ref_at_target_inject_location[valid_II]'*V

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

            # Normal output
            processed_S = target[region_S] .- refs_S_rotated*c

            # Model output:
            # processed_S .= targ_M .- vec(sum(reshape(c,1,:) .* refs_M_rotated,dims=2))

            phot_out = phot_of_target_at_inject_location .- dot(c, phot_S_rotated)
            
            throughput = phot_out/phot_of_target_at_inject_location

            # Throughput normalize (just once for full anulus)
            processed_S ./= throughput


            pixels = collect(StatsBase.trim(vec(processed_S); prop=0.05))
            if length(pixels) <= 0
                continue
            end
            clip_std = std(pixels) # Find std, while ignoring 10% most deviated pixels
            region_S_valid = all(
                (-4clip_std .< processed_S .< 4clip_std),
                dims=2
            )[:]
            # Don't do this clipping if we end up with nothing left, or take too much
            if count(region_S_valid) == 0 || mean(region_S_valid) < 0.25
                region_S_valid .= true
            end

            new_SNR = 1.0 ./ sqrt(mean(processed_S[region_S_valid].^2))

            print("\trotthreshpx=$rotthreshpx_val  \tn_SVD=$n_SVD \t: $(new_SNR)")
            if new_SNR > best_SNR
                print("\t*")
                best_SNR = new_SNR
                out[region_S] .= processed_S
            end
            println()
        end
    end

    println("improvement = $(best_SNR/init_SNR)")
    println()
end

function loci2_frame(fname, refnames_pattern, rotthreshpx, region_S, region_O; kwargs...)
    target = load(fname)
    outfname = replace(fname, ".fits"=>".loci.fits")
    if isfile(outfname)
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
        loci2_region!(out, target, refs, rotthreshpx, reg_S.>0, reg_O.>0, refcube; kwargs...)
    end
    AstroImages.writefits(outfname, out)
    return out
end



function loci2_all(fnames_pattern, rotthreshpx, regions_S, regions_O; force=false, kwargs...)
    fnames = glob(fnames_pattern)
    imgs = load.(fnames)
    if length(unique(size.(imgs))) > 1
        @error "Mismatched image dimensions"
        for (fname, img) in zip(fnames, imgs)
            println(fname, "\t", size(img)) 
        end
        error()
    end
    refcube=stack(imgs)
    i = 0
    for (fname, targ) in zip(fnames, imgs)
        i+=1
        ii = setdiff(axes(refcube,3),i)
        refs_this = @view imgs[ii]
        refcube_this = @view refcube[:,:,ii]

        outfname = replace(fname, ".fits"=>".loci.fits")
        if !force && isfile(outfname) &&  Base.Filesystem.mtime(outfname) > Base.Filesystem.mtime(fname)
            continue
            # out = load(outfname)
        else
            out = deepcopy(targ)
            fill!(out, NaN)
        end

        for (reg_S, reg_O) in zip(regions_S, regions_O)
            loci2_region!(out, imgs[i], refs_this, rotthreshpx, reg_S.>0, reg_O.>0, refcube_this; kwargs...)
        end
        AstroImages.writefits(outfname, out)      
        println(outfname, "\t($i)")  
    end
end