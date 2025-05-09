using AstroImages
using Glob
using CairoMakie:Makie

export seq2gif

function seq2gif(
    pattern::AbstractString,
    outname=replace(pattern, "*"=>"_", ".fits"=>".mp4",".gz"=>""); kwargs...
    )
    fnames = Glob.glob(pattern)
    return seq2gif(fnames, outname; kwargs...)
end
function seq2gif(
    fnames::AbstractArray{<:AbstractString},
    outname;
    clims=Percent(99.5),
    stretch=identity,
    crop=nothing,
    clims_per_frame=false
    )
    @info "loading frames"
    if isnothing(crop)
        imgs = stack(load.(fnames))
    else
        imgs = stack(map(fnames) do fname
            img = load(fname)
            if haskey(img, "STAR-X")
                cx = round(Int,img["STAR-X"])
                cy = round(Int,img["STAR-Y"])
            else
                cx = round(Int,mean(axes(img,1)))
                cy = round(Int,mean(axes(img,2)))
            end
            return img[cx-crop:cx+crop,cy-crop:cy+crop]
        end)
    end
    @info "rendering"
    if !clims_per_frame
        rendered = imview(imgs; clims, stretch)
    else
        @info "calculating clims per frame"
    end
    sz = size(imgs)[1:2]
    if min(sz...) < 400
        sz = (400,400)
    end
    fig = Makie.Figure(
        size=sz .* (1,1.05)
    )
    title = Makie.Observable(first(fnames))
    if !clims_per_frame
        frame = Makie.Observable(Makie.rotr90(@view rendered[:,:,1]))
    else
        frame = Makie.Observable(Makie.rotr90(imview(@view imgs[:,:,1];clims,stretch)))
    end
    ax = Makie.Axis(
        fig[1,1];
        title,
        aspect=Makie.DataAspect()
    )
    Makie.image!(ax,frame)

    Makie.record(fig, outname, axes(imgs,3); framerate=10) do i
        println("<- ", fnames[i])
        title[] = fnames[i]
        if !clims_per_frame
            frame[] = Makie.rotr90(@view rendered[:,:,i])
        else
            frame[] = Makie.rotr90(imview(@view imgs[:,:,i];clims,stretch))
        end
    end
    println(outname)
    return outname
end