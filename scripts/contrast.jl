using Statistics
using CairoMakie: Makie
using Glob

function imgsep(img::AstroImage)
    if haskey(img, "STAR-X")
        cx = img["STAR-X"]::Float64
        cy = img["STAR-Y"]::Float64
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
        cx = img["STAR-X"]::Float64
        cy = img["STAR-Y"]::Float64
    else
        cx = convert(Float64, mean(axes(img,1)))
        cy = convert(Float64, mean(axes(img,2)))
    end
    xs = axes(img,1) .- cx
    ys = axes(img,2) .- cy
    θ = atan.(ys, xs)
    return θ
end
function contrast(img::AstroImage; step=4)
    rs = imgsep(img)
    stop = maximum(rs)
    bins = range(;start=step/2, step, stop)
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
            cont[i] = std(view(img, mask))
        end
    end
    return bins, cont
end
function contrast(pattern::AbstractString; force=false)
    fnames = Glob.glob(pattern)
    return contrast(fnames; force)
end
function contrast(fnames::AbstractVector{<:AbstractString}; force=false)
    contrasts = []
    for fname in fnames
        # Calculate contrast if not there, or load
        outfname_txt = replace(fname, "*"=>"_", ".fits"=>".contrast.txt", ".gz"=>"")
        if !force && isfile(outfname_txt)
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
        seps,conts = contrast(img)
        push!(contrasts, (seps, conts))
        open(outfname_txt, write=true) do f
            println(f, "sep_px, 1sigcont")
            for (sep,cont) in zip(seps,conts)
                println(f, "$sep, $cont")
            end
        end
        println(outfname_txt)
    end
    return contrasts
end
function contrastplot(pattern::AbstractString; force=false)
    fnames = Glob.glob(pattern)
    contrasts = contrast(fnames; force)
    outfname_img = replace(pattern, "*"=>"_", ".fits"=>".png", ".gz"=>"")
    fig = Makie.Figure(
        size=(900,600)
    )
    # TODO: could add separation in MAS as a second x axis.
    ax = Makie.Axis(
        fig[1,1],
        xlabel="separation [px]",
        ylabel="1σ contrast",
        yscale=Makie.pseudolog10
    )
    for (fname, (seps,conts)) in zip(fnames,contrasts)
        Makie.lines!(
            ax,
            seps,
            conts,
            label=fname
        )
    end
    Makie.Legend(fig[1,2], ax, position=:rt, fontsize=8)
    Makie.save(outfname_img, fig)
    fig
end