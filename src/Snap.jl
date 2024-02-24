raw"""
```
.     _____ _   _          _____                          .
.    / ____| \ | |   /\   |  __ \      .-.,="``"=.        .
.   | (___ |  \| |  /  \  | |__) |     '=/_       \       .
.    \___ \| . ` | / /\ \ |  ___/       |  '=._    |      .
.    ____) | |\  |/ ____ \| |            \     `=./`,     .
.   |_____/|_| \_/_/    \_\_|             '=.__.=' `='    .
```

This package is a pipeline and collection of tools for
processing direct-imaging sequences.
```
"""
module Snap

function banner()
    return  raw""".     _____ _   _          _____                          .
                  .    / ____| \ | |   /\   |  __ \      .-.,="``"=.        .
                  .   | (___ |  \| |  /  \  | |__) |     '=/_       \       .
                  .    \___ \| . ` | / /\ \ |  ___/       |  '=._    |      .
                  .    ____) | |\  |/ ____ \| |            \     `=./`,     .
                  .   |_____/|_| \_/_/    \_\_|             '=.__.=' `='    ."""
end

using LibGit2: LibGit2 # Used for saving the pipeline commit hash in image headers
using TOML:TOML
using Dates
using Printf
using AstroLib: AstroLib
using Statistics, StatsBase
using AstroImages
using StaticArrays
using Optim

include("config.jl")
include("header-normalize.jl")
include("bad-pixels.jl")
include("find-star.jl")
include("dewarp-nirc2.jl")
include("crop.jl")

const as_per_rad = 206265




end