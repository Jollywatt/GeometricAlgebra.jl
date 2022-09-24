#= Run this script interactively: `julia --project=docs/ -i make.jl`
... or with arguments `julia --project=docs/ make.jl [test|fix|make|deploy]` =#

cd(joinpath(".", dirname(@__FILE__)))
using Pkg; Pkg.activate(".")

using Documenter, Revise, GeometricAlgebra

const project_root = pathof(GeometricAlgebra) |> dirname |> dirname

# apply setup code to all doctests in doc strings
DocMeta.setdocmeta!(GeometricAlgebra, :DocTestSetup, quote
    using Revise, GeometricAlgebra
end; recursive=true)

make() = makedocs(
    sitename="GeometricAlgebra.jl",
    root=joinpath(project_root, "docs"),
    modules=[GeometricAlgebra],
    pages=[
        "index.md",
        "design.md",
        "reference.md",
    ],
)

deploy() = deploydocs(
    repo = "github.com/Jollywatt/GeometricAlgebra.jl.git",
)

test() = doctest(GeometricAlgebra)

fix() = begin
    Revise.revise()
    doctest(GeometricAlgebra, fix=true)
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
