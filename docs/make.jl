#= Run this script interactively: `julia --project=docs/ -i make.jl`
... or with arguments `julia --project=docs/ make.jl [test|fix|make|deploy]` =#

cd(joinpath(".", dirname(@__FILE__)))
using Pkg; Pkg.activate(".")

using Documenter, Revise, Multivectors

const project_root = pathof(Multivectors) |> dirname |> dirname

# apply setup code to all doctests in doc strings
DocMeta.setdocmeta!(Multivectors, :DocTestSetup, quote
    using Revise, Multivectors
end; recursive=true)

make() = makedocs(
    sitename="Multivectors.jl",
    root=joinpath(project_root, "docs"),
    modules=[Multivectors],
    pages=["index.md"],
)

deploy() = deploydocs(
    repo = "github.com/Jollywatt/Multivectors.jl.git",
)

test() = doctest(Multivectors)

fix() = begin
    Revise.revise()
    doctest(Multivectors, fix=true)
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
