function prepregions(fname::AbstractString; kwargs...)
    img = load(fname)
    return prepregions(img; kwargs...)
end

"""
Prepare FITS files with two sets of region masks:
* subtraction regions
* optimization regions

Each HDU of the files will be another region.
"""
function prepregions(
    img::AstroImage;
    sub_inner_px=20,
    sub_thick_px=20,
    sub_outer_px=100,
    opt_inner_px=20,
    opt_thick_px=40,
    buf_px=10,
    num_sectors=3
)

    cx = convert(Float64, img["STAR-X"])
    cy = convert(Float64, img["STAR-Y"])
    xs = axes(img,1) .- cx
    ys = axes(img,2) .- cy
    Xs = collect(xs .+ 0 .* ys')
    Ys = collect(0 .* xs .+ ys')
    rs = sqrt.(xs.^2 .+ ys'.^2)
    θs = rem2pi.(atan.(ys', xs), RoundDown)

    # Don't include a new subtraction region if there is only 20% or less of a region left before
    # the outer separation.
    annuli_starts = range(start=sub_inner_px, stop=sub_outer_px-sub_thick_px/5, step=sub_thick_px)
    annuli_stops = range(start=sub_inner_px+sub_thick_px, stop=sub_outer_px, step=sub_thick_px)
    if length(annuli_stops ) < length(annuli_starts)
        annuli_starts = annuli_starts[begin:end-1]
    end

    sector_widths = 2π/num_sectors
    sector_starts = range(0, stop=2π-sector_widths, length=num_sectors)
    sector_stops = sector_starts .+ sector_widths
    
    S_masks = []
    O_masks = []
    for i_radius in eachindex(annuli_starts)
        for i_sector in 1:num_sectors
            S_mask = (
                (annuli_starts[i_radius] .<= rs .<= annuli_stops[i_radius]) .& 
                (sector_starts[i_sector] .<= θs .<= sector_stops[i_sector])
            )
            push!(S_masks, UInt8.(S_mask))
            # Optiization mask is everything within `opt_thick_px` (while greater than
            # `opt_inner_px`) or at same separation.
            # TODO: math needs double checking
            O_mask = (
                (opt_inner_px .< rs) .&  (
                    (
                        (annuli_starts[i_radius] .- opt_thick_px .< rs .<= annuli_starts[i_radius] .- buf_px) .|
                        (annuli_stops[i_radius] .+ buf_px .< rs .<= annuli_stops[i_radius] .+ opt_thick_px)
                    ) .|
                    (
                        (annuli_starts[i_radius] .- buf_px .<= rs .<= annuli_stops[i_radius] .+ buf_px) .& 
                        .!(sector_starts[i_sector] .<= θs .<= sector_stops[i_sector])
                    )
                )
            )
            push!(O_masks, UInt8.(O_mask))

            # O_mask = annuli_starts[i_radius] .<= rs .<= annuli_stops[i_radius]
            # # Find pixels close enough
            # for i_px in LinearIndices(img)
            #     for j_px in LinearIndices(img)
            #         dist = sqrt((Xs[i_px]-Xs[j_px])^2 + (Ys[i_px]-Ys[j_px])^2) 
            #         if dist < opt_thick_px
            #             O_mask[i_px] = true
            #         end
            #     end
            # end
            # O_mask .&= rs .> opt_inner_px 
            # # Omit pixels too close
            # for i_px in eachindex(img)
            #     for j_px in eachindex(img)
            #         dist = sqrt((Xs[i_px]-Xs[j_px])^2 + (Ys[i_px]-Ys[j_px])^2) 
            #         if dist < buf_px
            #             O_mask[i_px] = false
            #         end
            #     end
            # end
        end
    end
    AstroImages.writefits("masks.S.fits", S_masks...)
    AstroImages.writefits("masks.O.fits", O_masks...)
    return S_masks, O_masks
end