
export filetable

# Can be a list of paths, or a shell glob pattern.
function loadfiles(paths::AbstractArray)
    read.(paths, RasterImage)
end
function loadheaders(paths::AbstractArray)
    map(paths) do path
        FITS(path) do fit
            read_header(first(fit))
        end
    end
end
loadfiles(pathspec::AbstractString) = loadfiles(glob(pathspec))
loadheaders(pathspec::AbstractString) = loadfiles(glob(pathspec))
function loadheader(path::AbstractString, hdunum=1)
    FITS(path) do fit
        read_header(fit[hdunum])
    end
end


"""
    filetable("sequence/profile.toml", save=true, verbose=true)

This prints out a listing of the fits files in the `raw` directory
of the config file in a nice table format.

You can use it to identify the right frames, and just copy them
right into the .toml config file under the right section.

If `save=true`, then write the file table to a text file to 
`sequence/cal/profile.filelisting.txt`. If `verbose=true`, then
print it to the screen (both can be true).

Returns: nothing.
"""
function filetable(conf_fname; verbose=true, save=true, showplots=false)

    # Read the configuration
    conf = SNAP.readconfig(conf_fname)
    cal = conf["calibrate"]
    λoverD = conf["telescope"]["λoverD"];

    # Create a nice sanitized name for reporting these results
    if !occursin('/', conf_fname)
        conf_fname = joinpath(".", conf_fname)
    end
    profile = replace(conf_fname, ".toml" => "")
    profilename = splitpath(profile)[end-1]
    try
        profilename = replace(joinpath(splitpath(profile)[end-1:end]...), r"\W"=>"-")
    catch
    end
    profdir = dirname(profile)
    caldir = joinpath(profdir,"cal",basename(profile))
    if !isdir(joinpath(profdir,"cal"))
        verbose && @info("Creating `cal` directory next to profile", caldir)
        mkpath(joinpath(profdir,"cal"))
    end

    if haskey(conf["calibrate"], "base_path")
        prefix = conf["calibrate"]["base_path"]
    else
        prefix = joinpath(profdir, "raw")
    end
    fnames = [
        glob("*.fits*", prefix);
        glob("*/*.fits*", prefix);
    ]
    verbose && @info("Found $(length(fnames)) fits files")

    if save
        outfile = open("$caldir.filelisting.txt", "w")
        verbose && @info("Writing listing to $caldir.filelisting.txt")
    end

    table = (;
        paths=String[],
        date=String[],
        dims=[],
        totexp=[],
        object=String[],
        target=String[]
    )
    
    for fname in fnames
        try
            headers = loadheader(fname)

            # Keck NIRC2
            if haskey(headers, "TELESCOP") &&
                occursin("Keck", headers["TELESCOP"]) #&&
                # occursin("NIRC2", headers["CURRINST"])
                s = headers["NAXIS1"], headers["NAXIS2"]
                d = headers["DATE-OBS"]
                o = headers["OBJECT"]
                t = headers["TARGNAME"]
                t = headers["TARGNAME"]
                itime = headers["ITIME"]
                coadds = headers["COADDS"]
                totexp = itime*coadds
                wav = headers["CENWAVE"]
            # Gemini South or North (after GPI 2 upgrade)
            elseif haskey(headers, "TELESCOP") &&
                occursin("Gemini", headers["TELESCOP"])

                headers2 = loadheader(fname,2)

                s = headers2["NAXIS1"], headers2["NAXIS2"]
                d = headers["DATE-OBS"]*headers["UT"]
                o = headers["OBJECT"]
                t = headers["OBSTYPE"]
                totexp = headers2["ITIME"]
                wav = ""
            else
                @error("Unrecognized instrument. Assuming DATE-OBS, OBJECT, TARGNAME, and TOTEXP are present", maxlog=1)
                s = (
                    haskey(headers, "NAXIS1") ? headers["NAXIS1"] : "?",
                    haskey(headers, "NAXIS2") ? headers["NAXIS2"] : "?"
                )
                d = haskey(headers, "DATE-OBS") ? headers["DATE-OBS"] : "?"
                o = haskey(headers, "OBJECT") ? headers["OBJECT"] : "?"
                t = haskey(headers, "TARGNAME") ? headers["TARGNAME"] : "?"
                totexp = haskey(headers, "TOTEXP") ? headers["TOTEXP"] : "?"
                wav = ""
            end
            fnames = splitpath(fname)
            if fnames[1] == "sequences"
                popfirst!(fnames)
                popfirst!(fnames)
            end
            # Note: not using joinpath here on purpose, standardize output with forward slash.
            fname = join(fnames, "/")

            push!(table.paths, fname)
            push!(table.date, d)
            push!(table.dims, s)
            push!(table.totexp, totexp)
            push!(table.object, o)
            push!(table.target, t)

            verbose && println("\"$fname\", #\t$d\t$s\t$itime\t$coadds\t$o\t$t\t$wav")
            save && println(outfile, "\"$fname\", #\t$d\t$s\t$itime\t$coadds\t$o\t$t\t$wav")
        catch err
            # rethrow(err)
            println("#", typeof(err))
        end
    end
    save && close(outfile)

    return table
end

