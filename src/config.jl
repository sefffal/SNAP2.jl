
"""
    readconfig(filepath::AbstractString)
"""
function readconfig(filepath::T) where {T <: AbstractString}

    config = TOML.parsefile(filepath)
    config["telescope"]["λoverD"] =
        config["telescope"]["λ"] / config["telescope"]["D"] *
        as_per_rad * 1000
    return config

end