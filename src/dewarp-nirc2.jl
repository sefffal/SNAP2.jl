using CoordinateTransformations: Transformation
using ImageTransformations: ImageTransformations

const nirc2_distort_X_post20150413_v1 = load(joinpath(@__DIR__, "data", "nirc2_distort_X_post20150413_v1.fits"))
const nirc2_distort_Y_post20150413_v1 = load(joinpath(@__DIR__, "data", "nirc2_distort_Y_post20150413_v1.fits"))
const nirc2_distort_X_pre20150413_v2 = load(joinpath(@__DIR__, "data", "nirc2_distort_X_pre20150413_v2.fits"))
const nirc2_distort_Y_pre20150413_v2 = load(joinpath(@__DIR__, "data", "nirc2_distort_Y_pre20150413_v2.fits"))

function nirc2_dewarp!(img::AstroImage)
    # TODO: adapt
    after_20150413 = true
    try
        date_obs = parse(Date, img["DATE-OBS"])
        after_20150413 = date_obs > parse(Date, "2015-04-13")  
    catch
        @warn "Could not parse DATE-OBS header, assuming after 2015-04-13 for distortion correction" maxlog=5
        after_20150413 = true
    end
    if after_20150413
        dx = nirc2_distort_X_post20150413_v1
        dy = nirc2_distort_Y_post20150413_v1
    else
        dx = nirc2_distort_X_pre20150413_v2
        dy = nirc2_distort_Y_pre20150413_v2
    end

    dx = crop_A_to_B(dx, img)
    dy = crop_A_to_B(dy, img)
        
    tfrm = DistortionTransformation(dx, dy)

    img .= ImageTransformations.warp(img, tfrm, axes(img))
    img["DISTCORR"] =  (after_20150413 ? "post20150413_v1" : "pre20150413_v2")
    push!(img, History, "$(Date(Dates.now())): Applied distortion correction")

    return img

end



"""

Transformation object for applying an distortion correction map
to e.g. an image.
"""
struct DistortionTransformation{T} <: Transformation where {T <: AbstractMatrix} 
    dx::T
    dy::T
end
struct InvDistortionTransformation{T} <: Transformation where {T <: AbstractMatrix} 
    dx::T
    dy::T
end


(tfrm::DistortionTransformation)(x,y) = tfrm(SVector(x,y))
function (tfrm::DistortionTransformation)(position)
# function (tfrm::DistortionTransformation)(position::AbstractVector{<:Integer})

    # @info "type" typeof(position) position maxlog=20
    # Given: pixel coordinates `position`
    # Need: pixel coordinates after shifting due to distortion map.

    # Question: can we trust the positions will always be integers?
    # Start with yes, then consider adding interpolation into the distortion maps.
    i,j = position
    i = min(max(i, firstindex(tfrm.dx,1)), lastindex(tfrm.dx,1))
    j = min(max(j, firstindex(tfrm.dx,2)), lastindex(tfrm.dx,2))
    δx = tfrm.dx[i,j]
    δy = tfrm.dy[i,j]

    return position[1:2] .+ SVector(δx, δy)
end

(tfrm::InvDistortionTransformation)(x,y) = tfrm(SVector(x,y))
function (tfrm::InvDistortionTransformation)(position)
# function (tfrm::DistortionTransformation)(position::AbstractVector{<:Integer})

    # @info "type" typeof(position) position maxlog=20
    # Given: pixel coordinates `position`
    # Need: pixel coordinates after shifting due to distortion map.

    # Question: can we trust the positions will always be integers?
    # Start with yes, then consider adding interpolation into the distortion maps.
    i,j = position
    i = min(max(i, firstindex(tfrm.dx,1)), lastindex(tfrm.dx,1))
    j = min(max(j, firstindex(tfrm.dx,2)), lastindex(tfrm.dx,2))
    δx = tfrm.dx[i,j]
    δy = tfrm.dy[i,j]

    return position .- SVector(δx, δy)
end

import Base.inv
function Base.inv(tfrm::DistortionTransformation)

    # Given: DistortionTransformation
    # Need: Inverse of DistortionTransformation, i.e. the transform that would undo it.

    # Niave approach: take negative of distortion maps.
    # Downsides: works for small distortions, fails for big ones >0.5 pix.

    # Try to find a better approach
    # If I'm at pixel x=10, and the distortion map says pixel 8 -> pixel 10,
    # then I need to go back to 8. Note that pixel 9 could be different right?
    # Maybe we make that assumption...
    # There is no general way to handle this, so follow the niave approach.

    # Looking at the NIRC2 distortion maps, the gradient in the distortion is at most 0.3
    # so this is actually not a bad assumption.
    inv_tfrm = InvDistortionTransformation(tfrm.dx, tfrm.dy)

    return inv_tfrm
end