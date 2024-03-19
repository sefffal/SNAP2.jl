export snrmap
# function snrmap(pattern::AbstractString; force=false,step=4)
#     fnames = Glob.glob(pattern)
#     for fname in fnames
#         outfname = replace(fname,".fits"=>".snr.fits")
#         if !force && isfile(outfname) && Base.Filesystem.mtime(outfname) > Base.Filesystem.mtime(fname)
#             continue
#         end
#         img = load(fname)
#         # Get contrast curve
#         sep, cont = only(contrast(fname; force,step))
#         contrast_interpolator = Interpolations.linear_interpolation(
#             sep, cont,
#             extrapolation_bc=Line()
#         )
#         r = imgsep(img)
#         snr = img ./ contrast_interpolator.(r)

#         push!(snr, History, "$(Date(Dates.now())): converted to SNR map (used cont curve).")
#         AstroImages.writefits(outfname, snr)
#         println(outfname, "\t snr")
#     end
# end

function snrmap(pattern::AbstractString; force=false,step=4)
    fnames = Glob.glob(pattern)
    for fname in fnames
        outfname = replace(fname,".fits"=>".snr.fits")
        if !force && isfile(outfname) && Base.Filesystem.mtime(outfname) > Base.Filesystem.mtime(fname)
            continue
        end
        img = load(fname)
        # Get contrast curve
        contmap = only(contrastmap(fname; force,step))
        snr = img ./ contmap

        push!(snr, History, "$(Date(Dates.now())): converted to SNR map (used cont map).")
        AstroImages.writefits(outfname, snr)
        println(outfname, "\t snr")
    end
end