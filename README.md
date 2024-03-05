# SNAP2
```
.     _____ _   _          _____                          .
.    / ____| \ | |   /\   |  __ \      .-.,="``"=.        .
.   | (___ |  \| |  /  \  | |__) |     '=/_  v2   \       .
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

## General Process

The SNAP2 pipeline uses a template `.toml` file to select files for calibration. After that, each function takes file paths or shell "glob" patterns, eg. "n*.fits". 
Most steps only operate on new or changed files. This allows them to be used for live-processing.
If you want to re-run a step from scratch, pass `force=true`.

All steps generate temporary FITS files, so make sure that you have plently of disk space before starting!

After calibration, each step appends a new suffix to the end of the filename. For example,
`myseqname.n0101.cal.bgsub.medsub.rotnorth.fits` has been calibrated, background-subtracted, median-subtracted, and rotated North-up.

Operations that reduce several files take a "glob" pattern, eg. `myseqname.n*.cal.fits`. If the operation produces a single output file, eg due to stacking, the output filename will be the "glob" pattern with any `*` replaced by `_`, eg. `myseqname.n_.cal.stacked.fits`.

## Installation
1. [Install Julia 1.10+ via `juliaup`](https://julialang.org/downloads/): `curl -fsSL https://install.julialang.org | sh` (Mac and Linux).

2. Install SNAP. Open julia, and type  the character **`]`**. Then paste:
```
add https://github.com/sefffal/Snap2.jl.git
```

3. Start julia with multiple threads `julia --threads=auto`, and load SNAP: `using Snap`.

##  Pre-Processing

1. Create a template `.toml` file that assigns frames to use for eg. darks, flats, etc. There are also a few parameters to control pre-processing steps like background subtraction. **An example is listed at the bottom of this page.**

Starting with an blank `template.toml` file, run this function to print out a list of all your raw files. You can scroll through and copy-paste these entries into the right section of the `template.toml` file.
```julia
filetable("myseqname.toml")
```

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

## Quick median subtraction
Use this for eg live processing during an observing run. This function simply stacks the PSFs with a median, subtracts that median from all images, rotates them North-up, and stacks the residuals.
```julia
medsub("cal/abaur.n*.cal.bgsub.reg.fluxnorm.fits.gz")
```


## LOCI/KLIP subtraction

SNAP2 supports a LOCI / KLIP reduction. Note that these two algorithms are actually nearly indentical, despite the way they are described in the literature (Some implementations are referred to as LOCI, while others are referred to as KLIP). The only differ in the use of SVD versus PCA, which under some assumptions produce the same result.

Our implementation uses the SVD.

To perform LOCI PSF subtraction:

1. Generate subtraction regions. These are typically annular segments. The defaults below produce pairs of subtraction and optimization regions, separated by a buffer zone. This prevents over-subtraction of the planet PSF. Although I do not advise it, you can use a single subtraction zone by passing the same region list as both subtraction and optimization region arguments.

```julia
S_regions, O_regions = prepregions(
    "cal/myseqname.n0180.cal.bgsub.reg.fluxnorm.fits.gz",
    # Control the geometry of the subtraction and optimization regions.
    sub_inner_px=8,
    sub_outer_px=110,
    sub_thick_px=10,
    opt_inner_px=9,
    opt_thick_px=25,
    buf_px=4
);
```

2. Perform a LOCI subtraction on all frames. See `loci2_frame` and `loci2_region!` to only perform LOCI on a single frame, or single region of a frame (faster for experimentation)
```julia
loci2_all(
    "cal/myseqname.*.cal.bgsub.reg.fluxnorm.fits.gz",
    # Rejection distances in px to test for best SNR:
    0 : 2.5 : 10,
    S_regions,
    O_regions,
    # Full-width-at-half max of the planet PSF in pixels for SNR modelling.
    # Should be approximately 1.22 lambda/D.
    psf_fwhm=4.6,
    # Numpber of SVD components to use. 0 means all (no SVD truncation).
    N_SVD=[5,15,30,60,0]
)
# Optional: generate a movie
seq2gif("cal/myseqname.*.cal.bgub.reg.fluxnorm.loci.fits.gz")
```

3. Rotate subtracted images North-up and stack.
```julia
rotnorth("cal/myseqname.*.cal.bgsub.reg.fluxnorm.loci.fits.gz")
stackframes(median, "cal/myseqname.*.cal.bgsub.reg.fluxnorm.loci.rotnorth.fits.gz")
```

You could also perform a contrast-weighted stack like so:
```julia
stackframes_contweight(median, "cal/myseqname.*.cal.bgsub.reg.fluxnorm.loci.rotnorth.fits.gz")
```


## Multi-Target SNR Optimization


SNR Optimization is a generalization of LOCI/KLIP to account for planet PSF self-subtraction. Rather than rejecting frames that lead to too much planet subtraction, we model the flux contamination of each reference frame and find the linear combination of all frames (including overlapping frames) that maximizes the planet's SNR. **Note: this is not the same as optimizing the parameters of a LOCI reduction!** This operates one level lower in a single step.

In multi-target SNR optimization, we furthermore simultaneously reduce groups of about 10 contiguous frames. The optimizer works directly in the North-up rotated reference frame to optimize the **STACKED** non-linear SNR of that region. 

This is somewhat expensive, since we have to create a set of rotated reference images for each North-up target image, resulting in very large systems of equations. 
To reduce the scale of this problem somewhat, we apply SVD truncation to any reference images that have little or no flux contaminaion with the target images (less than 10%). For images with higher flux contamination, we bin together images that are more than 0.99 correlated and have less than 2% difference in the forward modelled flux contamination.

This approach could be considered a hybrid SVD/KIP and SNR-Opt.

Intructions:
1. Prepare subtraction and optimization regions
```julia
S_regions, O_regions = prepregions(
    "cal/myseqname.n*.cal.bgsub.reg.fluxnorm.inject.fits.gz",
    sub_inner_px=6,
    sub_thick_px=10,
    sub_outer_px=56,
    opt_inner_px=6,
    opt_thick_px=15,
    buf_px=5,
    num_sectors=3
);
```

2. Run mult-target SNR optimization. You can provide a list of different `N_SVD` values to perform the reduction with different numbers of SVD components. Fewer components will extrapolate better from the optimization to the subtraction region (less overfitting) but might not model the noise as well. A value of `0` means to include all reference images (no SVD truncation).
```julia
snropt_multitarg(
    "cal/myseqname.n*.cal.bgsub.reg.fluxnorm.fits.gz",
    S_regions,
    O_regions,
    # Full-width-at-half max of the planet PSF in pixels for SNR modelling.
    # Should be approximately 1.22 lambda/D in detector pixels.
    psf_fwhm=9.28,
    # Numpber of SVD components to use. 0 means all (no SVD truncation).
    # SVD is only applied to reference frames with less than 10% flux contamination.
    N_SVD=[5,15,30,60,120,0]
);
```

3. Stack frames. You don't have to rotate them North-up first, since that happens during the multi-target SNR optimization.
```julia
stackframes(median, "cal/myseqname.n*.cal.bgsub.reg.fluxnorm.rotnorth.snropt.fits.gz")
```


## Single-Target SNR Optimization

These functions perform SNR-optimization on a single target image which has not yet been rotated North-up (its basically just LOCI with the math changed out).

1. Prepare regions (see above)
```julia
S_regions, O_regions = prepregions(
    "cal/myseqname.n*.cal.bgsub.reg.fluxnorm.inject.fits.gz",
    sub_inner_px=6,
    sub_thick_px=10,
    sub_outer_px=56,
    opt_inner_px=6,
    opt_thick_px=15,
    buf_px=5,
    num_sectors=3
);
```

2. Run single-target SNR-Optimization on all frames. See `snropt_frame` to only process a single frame and `snropt_region!` to only process a single region of a frame.
```julia
snropt_all(
    "cal/myseqname.n*.cal.bgsub.reg.fluxnorm.inject.fits.gz",
    S_regions,
    O_regions,
    psf_fwhm=9.28/2,
    N_SVD=[5,15,30,]
);
```

3. Rotate North-up and stack:
```julia
rotnorth("cal/myseqname.*.cal.bgsub.reg.fluxnorm.snropt.fits.gz")
stackframes(median, "cal/myseqname.*.cal.bgsub.reg.fluxnorm.snropt.rotnorth.fits.gz")
```



## Calibration template file
To perform the initial calibration, you must assign your raw frames to the appropriate sections in a template file like in this example.

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

