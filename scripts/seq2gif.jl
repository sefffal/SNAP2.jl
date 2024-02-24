using AstroImages
using GIFImages
using Glob
using CairoMakie:Makie
function seq2gif(
    pattern::AbstractString,
    outname=replace(pattern, "*"=>"_", ".fits"=>".mp4",".gz"=>"");
    clims=Percent(99.5),
    crop=nothing
    )
    fnames = Glob.glob(pattern)
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
    rendered = imview(imgs; clims)
    sz = size(imgs)[1:2]
    if min(sz...) < 400
        sz = (400,400)
    end
    fig = Makie.Figure(
        size=sz .* (1,1.05)
    )
    title = Makie.Observable(first(fnames))
    frame = Makie.Observable(Makie.rotr90(@view rendered[:,:,1]))
    ax = Makie.Axis(
        fig[1,1];
        title,
        aspect=Makie.DataAspect()
    )
    Makie.image!(ax,frame)

    Makie.record(fig, outname, axes(rendered,3); framerate=10) do i
        println("<- ", fnames[i])
        title[] = fnames[i]
        frame[] = Makie.rotr90(@view rendered[:,:,i])
    end
end