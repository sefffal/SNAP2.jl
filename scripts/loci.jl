using FITSIO
using StatsBase
using TSVD # truncated SVD
function loci_region(fname, refnames_pattern::AbstractString, rotthreshpx, region_S, region_O; kwargs...)
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
    out = loci_region!(out, target, refs, rotthreshpx, region_S, region_O; kwargs...)
    AstroImages.writefits(outfname, out)
end
using LinearSolve
function loci_region!(out, target, refs::AbstractVector{<:AstroImage}, rotthreshpx, region_S, region_O, refcube=stack(refs); N_SVD=60)
    rs = imgsep(target)
    sep = mean(rs[region_S])
    angle = deg2rad(target["ANGLE_MEAN"]::Float64)
    angles = deg2rad.(getindex.(refs, "ANGLE_MEAN")::Vector{Float64})
    valid_II = findall(sep .* rem2pi.(angle .- angles, RoundDown) .> rotthreshpx)

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
    valid .&= all(-3clip_std .< view(refs_O,valid,:) .< 3clip_std,dims=2) .& (-3clip_std .< view(targ_O,valid) .< 3clip_std)

    if count(valid) == 0
        out[region_S] .= target[region_S]
        return out
    end
    

    ## SVD
    # N_SVD = 60 # 30
    U, s, V = tsvd(Float32.([refs_O[valid,:];  refs_S]), N_SVD)
    refs_O_rotated = U[1:count(valid),:]
    refs_S_rotated = U[count(valid)+1:end,:]
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



    # out[region_S] .= target[region_S] .- vec(sum(reshape(c,1,:) .* refs_S,dims=2))
    out[region_S] .= target[region_S] .- vec(sum(reshape(c,1,:) .* refs_S_rotated,dims=2))

    return out
end

function loci_frame(fname, refnames_pattern, rotthreshpx, region_S, region_O; kwargs...)
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
        out = loci_region!(out, target, refs, rotthreshpx, reg_S.>0, reg_O.>0, refcube; kwargs...)
    end
    AstroImages.writefits(outfname, out)
    return out
end



function loci_all(fnames_pattern, rotthreshpx, regions_S, regions_O; force=false, kwargs...)
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
            loci_region!(out, imgs[i], refs_this, rotthreshpx, reg_S.>0, reg_O.>0, refcube_this; kwargs...)
        end
        AstroImages.writefits(outfname, out)      
        println(outfname, "\t($i)")  
    end
end