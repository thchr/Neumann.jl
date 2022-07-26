using Documenter
using Neumann

# make sure we actually do `using Neumann` before calling doctests (as described in
# https://juliadocs.github.io/Documenter.jl/stable/man/doctests/#Module-level-metadata)
DocMeta.setdocmeta!(Neumann, :DocTestSetup, :(using Neumann); recursive=true)

makedocs(
    modules = [Neumann],
    sitename = "Neumann.jl",
    authors = "Thomas Christensen <tchr@mit.edu> and contributors",
    repo = "https://github.com/thchr/Neumann.jl/blob/{commit}{path}#L{line}",
    format=Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://thchr.github.io/Neumann.jl"
    ),
    pages = [
        "API"                   => "api.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo   = "github.com/thchr/Neumann.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing,
    push_preview = true
)