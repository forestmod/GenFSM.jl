# To build the documentation:
#    - julia --project="." make.jl
#    - empty!(ARGS); include("make.jl")
# To build the documentation without running the tests:
#    - julia --project="." make.jl preview
#    - push!(ARGS,"preview"); include("make.jl")

using Documenter
using GenFSM

push!(LOAD_PATH,"../src/")
makedocs(sitename="GenFSM.jl Documentation",
         pages = [
            "Index" => "index.md",
            "GenFSM" => "GenFSM.md",
            "ScenarioLoader" => "ScenarioLoader.md",
            "Res" => "Res.md",
         ],
         format = Documenter.HTML(prettyurls = false),
         warnonly = true
)
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/forestmod/GenFSM.jl.git",
    devbranch = "main"
)
