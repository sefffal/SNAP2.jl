
"""
    excludebadpixels!(img::AstroImage, σ=2, bound=1)

Replace bad pixels in an image with NaN.

Use the fact that these images are significantly oversampled to infer that strong
deviations from the surrounding pixels invariably means a detector issue.
i.e. a lambda/D is greater than a few pixels, so large changes below that pixel scale
can't be "real".

Suggestion: if you are working with pixel values that can be less than 0, make sure
to use only positive lowmedthresh and highmedthresh, or you will get confusing (but
correct) results.
"""
function excludebadpixels!(img::AstroImage, σ=2, bound=1)
    for I in CartesianIndices(img.data)
        x, y = Tuple(I)
        x1 = x - 1
        x2 = x + 1
        y1 = y - 1
        y2 = y + 1
        # Catch bounds errors and skip this pixel
        try
            img.data[x1 - bound:x2 + bound, y1 - bound:y2 + bound]
        catch e
            continue
        end
        pixels = collect(
            px
            for px = img.data[x1 - bound:x2 + bound, y1 - bound:y2 + bound] if isfinite(px)
        )
        if length(pixels) > 0
            localthresh = median(pixels)
            localstd = σ * std(pixels)
            if img.data[I] < (localthresh - localstd) ||
               img.data[I] > (localthresh + localstd)
                img.data[I] = NaN
            end
        else
            img.data[I] = NaN
        end
    end
end

"""
    interpbadpixels(image::AstroImage)

Replace all of the bad pixels in the given image
with a bilinear interpolation.

Edge pixels will be replaced with the overall image median -- in future this could be smarter.

Flag bad pixels using excludebadpixels! before calling this function to remove them.
Returns the number of pixels that were replaced.
"""
function interpbadpixels!(img::AstroImage)
    pixels = collect(d for d = img if isfinite(d))
    if length(pixels) == 0
        throw(ArgumentError("Cannot interpolate all NaN image"))
    end
    global_med = median(pixels)
    mat = zeros(4,4)
    col = zeros(4)
    coeffs = zeros(4)
    for I in CartesianIndices(img)
        x, y = Tuple(I)
        x1 = x - 1
        x2 = x + 1
        y1 = y - 1
        y2 = y + 1
        # Catch bounds errors and ignore
        if !(x1 ∈ axes(img,1) && x2 ∈ axes(img,1) &&
             y1 ∈ axes(img,2) && y2 ∈ axes(img,2))
             continue
        end
        if !isfinite(img[I])
            mat .= @SArray[
                1 x1 y1 (x1 * y1)
                1 x1 y2 (x1 * y2)
                1 x2 y1 (x2 * y1)
                1 x2 y2 (x2 * y2)
            ]
            patch = view(img, x1:x2, y1:y2)
            if any(!isfinite, patch)
                f = filter(isfinite, patch)
                if isempty(f)
                    med = global_med
                else
                    med = median(f)
                end
            end
            # Get surrounding pixels if finite (fallback to image median otherwise)
            col .= @SArray[
                isfinite(img[x1, y1]) ? img[x1, y1] : med
                isfinite(img[x1, y2]) ? img[x1, y2] : med
                isfinite(img[x2, y1]) ? img[x2, y1] : med
                isfinite(img[x2, y2]) ? img[x2, y2] : med
            ]
            # Solve the system of equations without inverting
            coeffs .= mat \ col
            interp =
                coeffs[1] + coeffs[2] * x + coeffs[3] * y + coeffs[4] * x * y
            img[I] = interp
        end
    end
    return img
end