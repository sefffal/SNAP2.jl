
"""
    findstar(image::AstroImage, (xmin px, xmax px), (xmin px, xmax px))

Find the centre of a star in an image given a small `xrange` and `yrange` surrounding it (and only it).
Return a named tuple with x and y positions of the centre.
"""
function findstar(image::AstroImage, x::Tuple, y::Tuple; fwhm=5.0)

    cropped = @view image[x[1]:x[2], y[1]:y[2]]

    Imax = Tuple(argmax(cropped))
    
    # This gives center in interpolated frame, interpolated frame may not be centred
    # Calculate the center of the real frame using center of interpolated frame, start of interpolated frame
    centrex = (Imax[1] - 1)
    centrey = (Imax[2] - 1)

    return (;
        x = centrex + x[1],
        y = centrey + y[1],
        amplitude = maximum(cropped.data),
        background = median(cropped.data),
    )
end


gauss2D(A, offset, μx, μy, σx, σy, x, y) =
    A * exp(-(((x - μx)^2) / (2 * σx^2) + ((y - μy)^2) / (2 * σy^2))) + offset

function centrefit_fixedstd(full::AbstractArray,  xs::Tuple=(1,size(full,1)), ys::Tuple=(1,size(full,2)); std, neg=false) 
    measured = @view full[xs[1]:xs[2], ys[1]:ys[2]]
    ValType = eltype(measured)
    model(A, offset, μx, μy, x, y) = gauss2D(A, offset, μx, μy, std, std, x, y)

    II = CartesianIndices(measured)
    xmax, ymax = size(measured)
    function objective(params)
        A, offset, μx, μy = params
        result = 0.0
        if μx > xmax || μy > ymax
            return ValType(Inf)
        end
        try # sig might be infeasible so we have to handle this case
            for I in II
                x, y = Tuple(I)
                est = model(A, offset, μx, μy, x, y)
                meas = measured[I]
                resid = (meas - est)^2
                # Skip NaN pixels in image - they give no weight.
                if isfinite(resid)
                    result += resid
                end
            end
        catch e
            # If the gaussian throws an error.
            result = ValType(Inf)
        end
        return log(result)
    end

    pixels = filter(isfinite, measured)
    if length(pixels) == 0
        throw(ArgumentError("Cannot optimize against all NaN/Inf image."))
    end

    # Inital guess: put in coordsinates of maximum point, and its height.
    startI = argmax(measured)
    guess_offset = mean(pixels)
    if neg
        startI = argmin(measured)
        guess_offset = maximum(pixels)
    end
    guess = [
        # A, offset, μx, μy
        # maximum(pixels) - mean(pixels),
        measured[startI],
        guess_offset,# mean(pixels),
        startI[1],
        startI[2],
    ]
    # test for an example starting point
    result = optimize(
        objective,
        guess,
        Optim.LBFGS(),
        Optim.Options(
            f_tol=1e-8,
            g_tol=1e-8,
            show_trace=false,
            time_limit=10.0, # Soft upper limit in seconds before returning best guess.
            allow_f_increases=true
        ),
        autodiff=:forward,
    )
    if !Optim.converged(result)
        @error("Centre-fit did not converge")
        display(result)
        # Return our basic guess instead of some 
        # badly converged result.
        result.minimizer .= [
            # A, offset, μx, μy
            # maximum(pixels) - mean(pixels),
            measured[startI],
            mean(pixels),
            startI[1],
            startI[2],
        ]
    end

    return (;
        x = result.minimizer[3] + xs[1] - 1,
        y = result.minimizer[4] + ys[1] - 1,
        amplitude = result.minimizer[1],
        background = result.minimizer[2],
    )
end