
"""
Special processing specific to NIRC2 images. This is called automatically
when running read("*.fits", RasterImage) when TELESCOP starts with Keck, and
CURRINST="NIRC2".
"""
function header_normalize_nirc2!(img)

    if haskey(img, "PA") && isnan(img["PA"])
        img["PA"] = 0.0
    end

    # Check if angles alrady computed
    angles = header_get_pa_angles(img)
    if length(angles) == 0
        angle_mean, angles = calculate_north_angle(img)
        # img["ANGLE_MEAN"] = angle_mean
    else
        angle_mean = SNAP.meandegrees(angles)
    end
    angles .+= (angle_mean - SNAP.meandegrees(angles))
    if isnan(angle_mean)
        angle_mean = 0.0
    end
    push!(img, History, "$(Date(Dates.now())): Re-calculated ANGLE_MEAN")
    img["ANGLE_MEAN"] = angle_mean


    if haskey(img, "EFFWAVE")
        img["WAVE"] = img["EFFWAVE"]
    elseif haskey(img, "CENWAVE")
        img["WAVE"] = img["CENWAVE"]
    end

    # Normalize the object key
    object = lowercase(img["SHRNAME"]) == "closed" ? "dark" :
        replace(lowercase(img["OBJECT"]), " " => "")
    img["OBJECT"] = object

    # # Ensure we are working with MAS units instead of pixel offsets on the axes
    # if step(xrange) == 1.0
    #     xrange = xrange .* PlateScale
    # end
    # if step(yrange) == 1.0
    #     yrange = yrange .* PlateScale
    # end

    push!(img, History, "$(Date(Dates.now())): Normalized headers for SNAP pipeline.")


    return img
end



float_or_parse_hexages(num::Number) = num
float_or_parse_hexages(str::AbstractString) = AstroLib.ten(str)

const observatory_lat = +19.82525
const observatory_lon = -155.468889
"""
Calculate and return the starting angle, angular smear, and mean angle to north
of an image from NIRC2 narrow cam, given its headers.

The calculation is adapted from pyKLIP.
"""
function calculate_north_angle(headers)
    zp_offset = -0.262 # From Service et al 2016
    rotator_mode = headers["ROTMODE"]
    rotator_position = headers["ROTPOSN"] # Degrees
    instrument_angle = headers["INSTANGL"] # Degrees

    pa_deg = 0.0
    if rotator_mode == "vertical angle"
        parang = headers["PARANTEL"]
        pa_deg = parang + rotator_position - instrument_angle + zp_offset
    elseif rotator_mode == "position angle"
        pa_deg = rotator_position - instrument_angle + zp_offset
    elseif rotator_mode == "stationary"
        # TODO: Handle case where the instrument rotator is stationary.
        return NaN, [NaN]
    else
        throw(ArgumentError("Unknown rotator mode " * rotator_mode))
    end

    # Keck rotator bug.
            # parang = headers["PARANG"] # Degrees
    # This does not appear to actually work
    diff = 0.0
    # if haskey(headers, "ROTNORTH")
    #     pa_deg_idl = -headers["ROTNORTH"]
    #     # diff += (pa_deg_idl - pa_deg)
    #     # println(diff)
    #     @info "Using ROTNORTH header from IDL pipeline" maxlog=1
    # end

    # Now calculate smear

    # Get info for PA smearing calculation.
    # epochobj = headers["DATE-OBS"]
    # name = headers["TARGNAME"]
    expref = headers["ITIME"]
    coaddref = headers["COADDS"]
    sampref = headers["SAMPMODE"]
    msrref = headers["MULTISAM"]
    xdimref = headers["NAXIS1"]
    # ydimref = headers["NAXIS2"]
    dec = float_or_parse_hexages(headers["DEC"]) + headers["DECOFF"]

    if haskey(headers, "TOTEXP")
        totexp = headers["TOTEXP"]
    else
        # Calculate total time of exposure (integration + readout).
        if sampref == 2
            totexp = (expref + 0.18 * (xdimref / 1024.0)^2) * coaddref
        end
        if sampref == 3
            totexp =
                (expref + (msrref - 1) * 0.18 * (xdimref / 1024.0)^2) * coaddref
        end
    end
    # tinteg = totexp # [seconds]
    totexp = totexp / 3600.0 # [hours]

    # Get hour angle at start of exposure.
    tmpahinit = float_or_parse_hexages(headers["HA"]) # [deg]
    ahobs = 24.0 * tmpahinit / 360.0 # [hours]

    if totexp * 3600 > 1 # If greater than 1 second...
        # Estimate vertical position angle at each second of the exposure.
        vp = Float64[]
        vpref = 0.0
        for j in 0:(3600*totexp-1)
            ahtmp = ahobs + (j + 1.0 + 0.001) / 3600.0 # hours
            # TODO: par_angle and observatory_latitude
            push!(vp, par_angle(ahtmp, dec, observatory_lat))
            if j == 0
                vpref = vp[1]
            end
        end

        # Handle case where PA crosses 0 <--> 360.
        vp[vp.<0] .+= 360
        vp[vp.>360.0] .+= 360

        if vpref < 0
            vpref += 360
        end
        if vpref > 360
            vpref -= 360
        end

        # TODO: This is some crazy code... Should use use meandegrees()
        # that uses atan() to avoid this.
        # Check that images near PA=0 are handled correctly.
        if any(vp .> 350) && any(vp .< 10)
            vp[vp.>350] .-= 360
        end
        vpmean = mean(angle for angle = vp if isfinite(angle))

        if vpmean < 60 && vpref > 350
            vpmean += 360
            vp .+= 360
        end
        pa_deg_mean = pa_deg + (vpmean - vpref)

        # angle_start = pa_deg
        # angle_smear = vpmean - vpref
        angle_mean = pa_deg_mean
        return angle_mean+diff, vp.+diff
    else # Total exposure less than one second - treat as no smear
        return pa_deg+diff, [pa_deg+diff]
    end
end

"""
Compute the parallactic angle, given hour angle (HA [hours]),
declination (dec [deg]), and latitude (lat [deg]).  Returns
parallactic angle in [deg].
Source: pyKLIP
"""
function par_angle(HA, dec, lat)

    HA_rad = deg2rad(HA * 15.0) # [hours] -> [rad]
    dec_rad = deg2rad(dec)   # [deg] -> [rad]
    lat_rad = deg2rad(lat)   # [deg] -> [rad]

    parallang =
        -atan(
            -sin(HA_rad),  # [rad]
            cos(dec_rad) * tan(lat_rad) - sin(dec_rad) * cos(HA_rad),
        )

    return rad2deg(parallang) # [deg]
end


# const PlateScale = 9.971



# Adaptation of Stan Metchev's lost linearity routine from IDL (found in Wayback machine)
# Original paper: https://iopscience.iop.org/article/10.1088/0067-0049/181/1/62/pdf

#     pro linearize_nirc2,fitslist
# ;
# ; Procedure to linearize NIRC2 (treated as a single detector).
# ; Run on all images before running anything else.  
# ;
# ; fitslist	- list of FITS files in 1024x1024 format to linearize
# https://web.archive.org/web/20160316114416/http://www.astro.sunysb.edu/metchev/AO/linearize_nirc2.pro
function linearize_nirc2_data!!(img::AstroImage)
    coadds = img["COADDS"]
    coeff = (1.001,-6.9e-6,-0.70e-10)
    if eltype(img) <: Integer
        img = float.(img)
    end
    norm = @. coeff[1]+coeff[2]*img.data/coadds+coeff[3]*(img.data/coadds)^2
    img .= img ./ norm .* coadds # Note: unlike original, we do not normalize by the number of coadds
    return img
end