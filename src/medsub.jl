
export medsub

"""
Given a single reference, typically and average stack of all frames,
subtract it from the set of matching frames, rotate those residuals
North-up, and then stack.
"""
function medsub(pattern::AbstractString; force=false)
    fnames = Glob.glob(pattern)
    @info "stacking master"
    stackframes(median, pattern; force)
    reference_img_fname = replace(pattern, "*"=>"_", ".fits"=>".stacked.fits")

    @info "subtracting"
    ref_img = collect(load(reference_img_fname))
    for fname in fnames
        outname = replace(fname,".fits"=>".medsub.fits")
        if !force && isfile(outname) && Base.Filesystem.mtime(outname) > Base.Filesystem.mtime(fname)
            continue
        end
        img = load(fname)
        img .-= ref_img
        push!(img, History, "$(Date(Dates.now())): subtracted reference image.")
        img["REFSUB"] = reference_img_fname
        AstroImages.writefits(outname, img)
        println(outname)
    end

    @info "stacking residuals (not North up)"
    medsub_pattern = replace(pattern, ".fits"=>".medsub.fits")
    stackframes(median, medsub_pattern; force)

    @info "rotating residuals North up"
    rotnorth(medsub_pattern; force)

    @info "stacking north up residuals"
    medsub_rotnorth_pattern = replace(pattern, ".fits"=>".medsub.rotnorth.fits")
    stackframes(median, medsub_rotnorth_pattern; force)
    # stackframes_contweight(median, medsub_rotnorth_pattern; force)
end