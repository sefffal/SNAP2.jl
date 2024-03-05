# SNAP2
```
.     _____ _   _          _____                          .
.    / ____| \ | |   /\   |  __ \      .-.,="``"=.        .
.   | (___ |  \| |  /  \  | |__) |     '=/_       \       .
.    \___ \| . ` | / /\ \ |  ___/       |  '=._    |      .
.    ____) | |\  |/ ____ \| |            \     `=./`,     .
.   |_____/|_| \_/_/    \_\_|             '=.__.=' `='    .
```
The Signal to Noise Analysis Pipeline, the Next Generation! (tm)

This repository contains a re-write of the SNAP pipeline, designed to be more modular and inspectable.

Current support:
* NIRC2
* LOCI/KLIP/TLOCI
* SNR-Opt
* Multi-target SNR Opt

Unsupported (see SNAP 1):
* IRDIS
* SPHERE IFS
* GPI
* NIRI

## Instructions: General Process

The SNAP2 pipeline uses a template `.toml` file to select files for calibration. After that, each function takes file paths or shell "glob" patterns, eg. "n*.fits". 
Most steps only operate on new or changed files. This allows them to be used for live-processing.
If you want to re-run a step from scratch, pass `force=true`.

All steps generate temporary FITS files, so make sure that you have plently of disk space before starting!

After calibration, each step appends a new suffix to the end of the filename. For example,
`myseqname.n0101.cal.bgsub.medsub.rotnorth.fits` has been calibrated, background-subtracted, median-subtracted, and rotated North-up.

Operations that reduce several files take a "glob" pattern, eg. `myseqname.n*.cal.fits`. If the operation produces a single output file, eg due to stacking, the output filename will be the "glob" pattern with any `*` replaced by `_`, eg. `myseqname.n_.cal.stacked.fits`.

##  Instructions: Pre-Processing

1. Create a template `.toml` file with the following information. This specifies which frames to use for eg. darks, flats, etc. There are also a few parameters to control pre-processing steps like background subtraction. An example is listed at the bottom of this page.

2. Run calibration
```julia
calibrate_nirc2("myseqname.toml",) # see example at the bottom of the page.
# Optional: generate a movie
seq2gif("cal/myseqname.*.cal.fits.gz")
```

3. Find the rough star position (updates headers). The arguments are "file pattern", "PSF template path", and "max star distance from centre (px)" for the search.
```julia
# Optional: generate a cropped movie.
findstar_nirc2("cal/myseqname.*.cal.fits.gz", "cal/unsat.psf-cored.fits.gz", 140,) 
seq2gif("cal/myseqname.*.cal.fits.gz", crop=100) 
```

4. Run background subtraction. 
```julia
bgsub("./myseqname.toml",)
seq2gif("cal/myseqname.*.cal.bgub.fits.gz", clims=Percent(95)) # Optional: generate a cropped movie.
```

5. Registration. Re-align the images so that the star is centred in all frames. Uses a template PSF and cross-correlation. If you are registering a coronagraphic sequence, the template can have the middle removed so that it only searches the surrounding speckle pattern.
```julia
register("cal/myseqname.*.cal.bgsub.fits.gz", "cal/unsat.psf-cored.fits.gz",)
# Optional: generate a cropped movie
seq2gif("cal/myseqname.*.cal.bgsub.reg.fits.gz", crop=150, clims=px -> quantile(px, 0.95) .* (-0.5, 1)) 
```

6. Flux normalization. 
```julia
fluxnorm("cal/myseqname.*.cal.bgsub.reg.fits.gz")
# Optiona: generate a movie
seq2gif("cal/myseqname.*.cal.bgsub.reg.fluxnorm.fits.gz", crop=150, clims=px -> quantile(px, 0.999) .* (-0.5, 1)) 
```



See inline comments for a description of each parameter.
```toml
# General telescope/instrument parameters
[telescope]
    # What instrument was this? Used to select the loading & calibration engine.
    inst="NIRC2"
    # Effective wavelength of the observations, meters.
    "λ" = 4.6705e-6
    # Telescope diameter, meters
    D = 10.0

[calibrate]

    # All file paths are relative to this file.

    # Path to a bad-pixel file. It should be the same size as your data
    # crop or larger. White=bad pixel, black=good pixel.
    badpix = "../nirc2-upgrade2024-bpm.fits"
    # Paths to dark files.
    dark = [
        "../2024-01-23/raw/n0034.fits",
        "../2024-01-24/raw/n0078.fits",
        "../2024-01-24/raw/n0079.fits",
        # ...
    ]
    # Paths to dark files for your unsaturated PSF images.
    # If you are using non-coronagraphic data that is unsaturated (
    # and therefore don't need separate unsaturated images) you 
    # can just repeat the same files here as under `dark`
    dark_psf = [
        "../2024-01-23/raw/n0034.fits", 
        "../2024-01-24/raw/n0078.fits", 
        "../2024-01-24/raw/n0079.fits", 
        # ...
    ]
    # Images to use for flat fielding. Must having matching `flatoff`
    # entries.
    flaton = [
        "../2024-01-24/raw/n0017.fits", 
        "../2024-01-24/raw/n0018.fits", 
        "../2024-01-24/raw/n0019.fits",    
        # ..
    ]
    # Flat darks matching exposure and coadds of `flaton` entries.
    flatoff = [
        # ..
    ]
    # Pure sky images, eg where telescope has moved and the star is no
    # longer in the field of view.
    # These are used in bg-subtraction (in addition to any chopped images)
    sky = [
        "./raw/n0166.fits",
        "./raw/n0167.fits",
        "./raw/n0168.fits",
        # ...
    ]
    # Optional sky images to use for BG subtracting the unsaturated PSF
    # template images. These can usually be omitted.
    sky_psf = [
    ]
    # Path to unsaturated PSF template images, where star is not behind
    # a coronagraph. 
    # If you have a non-coronagraphic and unsaturated sequence, these
    # can just be some of the main science raw files (see below.)
    raw_psf_path = [
        "./raw/n0099.fits",
        "./raw/n0100.fits",
        # ...
    ]
    # Your actual data, where the star is in the frame.
    # List all files here.
    raw_path = [
        "./raw/n0101.fits",
        "./raw/n0102.fits",
        "./raw/n0103.fits",
        # ...
    ]

    # Where should we output the unsaturated template PSF?
    output_psf_path = "cal/unsat.psf.fits.gz"

    ## BG Subtraction

    # Mask star to this radius for background sub.
    # If chopping is used, only images more than half this separation
    # will be included in the least squares.
    sky_psf_radius_px = 50.0

    # Pick the best N other light images (far enough away) for chopping sky subtraction
    sky_N_most_similar = 100

    # Apply chopping based sky subtraction (vs only using sky frames listed above)
    chopped = true


    ## Finding the star

```



regions = prepregions("cal/myseqname.n0180.cal.bgsub.reg.fluxnorm.fits.gz", sub_thick_px=10, sub_outer_px=110, sub_inner_px=8, opt_inner_px=9, opt_thick_px=25, buf_px=4);
# loci_all("cal/myseqname.*.cal.bgsub.reg.fluxnorm.fits.gz", 7, regions[1], regions[2], force=true)
loci2_all("cal/myseqname.*.cal.bgsub.reg.fluxnorm.fits.gz", 5.6 * 0.65, regions[1], regions[2], force=true, psf_fwhm=4.6, N_SVD=30)
seq2gif("cal/myseqname.*.cal.bgub.reg.fluxnorm.sub.fits.gz")

rotnorth("cal/myseqname.*.cal.bgsub.reg.fluxnorm.sub.fits.gz", force=true)
stackframes(median, "cal/myseqname.*.cal.bgsub.reg.fluxnorm.sub.rotnorth.fits.gz", force=true)
