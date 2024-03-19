using Statistics
using CairoMakie: Makie
using Glob
using StatsBase

export imgsep, imgang, contrast, contrastmap, contrastplot

function imgsep(img::AstroImage)
    if haskey(img, "STAR-X")
        cx = Float64(img["STAR-X"])::Float64
        cy = Float64(img["STAR-Y"])::Float64
    else
        cx = convert(Float64, mean(axes(img,1)))
        cy = convert(Float64, mean(axes(img,2)))
    end
    xs = axes(img,1) .- cx
    ys = axes(img,2) .- cy
    rs = sqrt.(xs.^2 .+ ys'.^2)
    return rs
end
function imgang(img::AstroImage)
    if haskey(img, "STAR-X")
        cx = Float64(img["STAR-X"])::Float64
        cy = Float64(img["STAR-Y"])::Float64
    else
        cx = convert(Float64, mean(axes(img,1)))
        cy = convert(Float64, mean(axes(img,2)))
    end
    xs = axes(img,1) .- cx
    ys = axes(img,2) .- cy
    θ = atan.(ys', xs)
    return θ
end

#= 
Contrast Curve calculation
=#
function contrast(img::AstroImage; step=4)
    rs = imgsep(img)
    # stop = maximum(rs)
    maxfinite = maximum(rs .* isfinite.(img))
    bins = range(;start=step/2, step, stop=maxfinite)
    OutType= promote_type(eltype(img), Float16)
    cont = zeros(OutType, length(bins))
    mask = falses(size(img))
    mask_finite = isfinite.(img)
    for i in eachindex(bins)
        sep = bins[i]
        mask .= mask_finite .& (sep-step/2 .<= rs .<= sep+step/2)
        if count(mask) == 0
            cont[i] = NaN
        else
            # reject 5% of outlying pixels when calculating
            trimmed = StatsBase.trim(view(img, mask), prop=0.05)
            cont[i] = std(trimmed)
        end
    end
    return bins, cont
end
function contrast(pattern::AbstractString; force=false, step=4)
    fnames = Glob.glob(pattern)
    return contrast(fnames; force, step)
end
function contrast(fnames::AbstractVector{<:AbstractString}; force=false, step=4)
    contrasts = []
    for fname in fnames
        # Calculate contrast if not there, or load
        outfname_txt = replace(fname, "*"=>"_", ".fits"=>".contrast.txt", ".gz"=>"")#, ".rotnorth."=>".", ".rotback."=>".")
        # TODO: check itime
        if !force && isfile(outfname_txt) && Base.Filesystem.mtime(outfname_txt) > Base.Filesystem.mtime(fname)
            open(outfname_txt) do f
                header = readline(f) # skip first line
                seps = Float64[]
                conts = Float64[]
                for line in readlines(f)
                    sep,cont = split(line, ',')
                    push!(seps, parse(Float64, sep))
                    push!(conts, parse(Float64, cont))
                end
                push!(contrasts, (seps, conts))
            end
            continue
        end
        img = load(fname)
        seps,conts = contrast(img; step)
        push!(contrasts, (seps, conts))
        open(outfname_txt, write=true) do f
            println(f, "sep_px, 1sigcont")
            for (sep,cont) in zip(seps,conts)
                println(f, "$sep, $cont")
            end
        end
        println(outfname_txt)
    end
    return identity.(contrasts) # promote container type
end


function contrastplot(pattern::AbstractString; force=false, platescale=9.971, step=4)
    fnames = Glob.glob(pattern)
    outfname_img = replace(pattern, "*"=>"_", ".fits"=>".png", ".gz"=>"")
    fig = contrastplot(fnames, outfname_img; force, platescale, step)
    Makie.save(outfname_img, fig)
    return fig
end
function contrastplot(fnames::AbstractArray{<:AbstractString}; force=false, platescale=9.971, step=4)
    contrasts = contrast(fnames; force, step)
    fig = Makie.Figure(
        size= length(fnames) > 1 ? (700,800) : (700,600)
    )
    ax = Makie.Axis(
        fig[1,1],
        xlabel="separation [px]",
        ylabel="1σ contrast",
        yscale=log10,
        xgridvisible=false,
        yminorticksvisible=true,
        yminorticks=Makie.IntervalsBetween(9),
        yticks=10.0 .^ (-10:10),
        xminorticks=Makie.IntervalsBetween(10),
        xminorticksvisible=true,
        ygridvisible=true,
        yminorgridvisible=true,
    )
    maxsep = 0
    for (fname, (seps,conts)) in zip(fnames,contrasts)
        Makie.scatterlines!(
            ax,
            seps,
            conts,
            label=fname
        )
        maxsep = max(maxsep, maximum(seps))
    end
    Makie.xlims!(ax, 0, maxsep*1.02)
    ax2 = Makie.Axis(
        fig[1,1],
        xaxisposition = :top,
        # yaxisposition = :right,
        xlabel="separation [mas]",
        ylabel="1σ contrast",
        yscale=log10,
        xgridvisible=false,
        yminorticksvisible=true,
        yminorticks=Makie.IntervalsBetween(9),
        yticks=10.0 .^ (-10:10),
        xminorticks=Makie.IntervalsBetween(10),
        xminorticksvisible=true,
        ygridvisible=true,
        yminorgridvisible=true,
    )
    Makie.linkyaxes!(ax, ax2)
    Makie.xlims!(ax2, 0, maxsep*platescale)
    if length(fnames) > 1
        Makie.Legend(fig[2,1], ax, position=:rt, fontsize=8, tellheight=true, tellwidth=false)
    end
    fig
end



#= 
Contrast Map calculation
At each pixel mask out a circular region and compute contrast as normal.
The result is a map.
=#
function contrastmap(img::AstroImageMat; step=4, maskr=15)
    rs = imgsep(img)
    θs = imgang(img)
    xs = Matrix{Float64}(rs .* cos.(θs))
    ys = Matrix{Float64}(rs .* sin.(θs))

    mask = falses(size(img))
    mask_finite = isfinite.(img)
    out = similar(img)
    fill!(out, NaN)

    if all(!isfinite, mask_finite)
        @error "all non-finite image"
        return out
    end
    # i = 0
    # l = count(mask_finite)
    # Threads.@threads 
    for I in CartesianIndices(img)[mask_finite]#[begin:100:end]
        # i += 1
        xi, yi = I[1], I[2]
        x = xs[I]
        y = ys[I]
        sep = rs[xi,yi]
        mask .= mask_finite
        mask[mask] .&=  @views (
            sep-step/2 .<= rs[mask] .<= sep+step/2
        ) 
        if count(mask) == 0
            continue
        end
        mask[mask] .&= @views sqrt.((xs[mask] .- x).^2 + (ys[mask] .- y).^2) .> maskr
        # display(imview(mask))
        if count(mask) == 0
            continue
        end
        # reject 5% of outlying pixels when calculating
        trimmed = StatsBase.trim(view(img, mask), prop=0.05)
        out[I] = std(trimmed)
    end
    
    return out
end

function contrastmap(pattern::AbstractString; force=false, step=4)
    fnames = Glob.glob(pattern)
    return contrastmap(fnames; force, step)
end
function contrastmap(fnames::AbstractVector{<:AbstractString}; force=false, step=4)
    contrasts = []
    for fname in fnames
        # Calculate contrast if not there, or load
        outfname = replace(fname, ".fits"=>".contrast.fits")
        # TODO: check itime
        if !force && isfile(outfname) && Base.Filesystem.mtime(outfname) > Base.Filesystem.mtime(fname)
            push!(contrasts, load(outfname))
            continue
        end
        img = load(fname)
        cont = contrastmap(img; step)
        AstroImages.write(outfname, cont)
        push!(contrasts, cont)
        println(outfname)
    end
    return identity.(contrasts) # promote container type
end
