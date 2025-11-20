#= Run this script interactively: `julia -i make.jl`
... or with arguments `julia make.jl [test|fix|make|deploy]` =#

using Documenter, DocumenterCitations
using Revise
using LiveServer

using GeometricAlgebra
using GeometricAlgebra.MiniCAS

const project_root = pathof(GeometricAlgebra) |> dirname |> dirname

# apply setup code to all doctests in doc strings
DocMeta.setdocmeta!(GeometricAlgebra, :DocTestSetup, quote
    using Revise
    using GeometricAlgebra
    using GeometricAlgebra.MiniCAS
end; recursive=true)

const bibfile = joinpath(project_root, "docs/src/references.bib")

function make(; kwargs...)
    macros = Dict(
        "\\lcontr" => "\\rfloor",
        "\\rcontr" => "\\lfloor",
    )

    makedocs(
        sitename="GeometricAlgebra.jl",
        root=joinpath(project_root, "docs"),
        modules=[GeometricAlgebra],
        plugins=[CitationBibliography(bibfile)],
        pages=[
            "index.md",
            "design.md",
            "Mathematical Background" => [
                "theory/basics.md",
                "theory/products.md",
                "theory/automorphisms.md",
                "theory/dualities.md",
                "theory/rotors.md",
                "theory/special.md",
                "theory/references.md",
            ],
            "Models of geometry" => [
                "theory/cga.md",
            ],
            "docstrings.md",
        ],
        format=Documenter.HTML(
            mathengine=KaTeX(Dict(
                :macros => macros,
                :minRuleThickness => 0.06,
            )),
            assets = ["assets/favicon.ico"],
        ),
        remotes=Dict(".." => Remotes.GitHub("Jollywatt", "GeometricAlgebra.jl"));
        kwargs...
    )
end
draft() = make(; draft = true)

deploy() = deploydocs(
    repo = "github.com/Jollywatt/GeometricAlgebra.jl.git",
)

test() = doctest(GeometricAlgebra)

fix() = begin
    Revise.revise()
    doctest(GeometricAlgebra, fix=true, plugins=[CitationBibliography(bibfile)])
    nothing
end

if isempty(ARGS)
    @info """Run this script in interactive mode, and call:
         - `make()` to build documentation locally
         - `deploy()` to build and deploy to github
         - `test()` to run doctests
         - `fix()` to fix doctests
        Keep the session alive; changes will be revised and successive runs will be faster.
        Alternatively, run this script passing a command (without parentheses) as an argument.
        """
    if !isinteractive()
        make()
        deploy()
    end
end

if "test" in ARGS
    test()
end

if "fix" in ARGS
    fix()
end

if "deploy" in ARGS
    make()
    deploy()
elseif "make" in ARGS
    make()
end

function cleardocs!(mod, name)
    Docs.meta(mod)[Docs.Binding(mod, name)] = Docs.MultiDoc()
end

cd(joinpath(project_root, "docs")) do
    make()
end
