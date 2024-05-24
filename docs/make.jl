using Documenter
using GenFSM

push!(LOAD_PATH,"../src/")
makedocs(sitename="GenFSM.jl Documentation",
         pages = [
            "Index" => "index.md",
            "An other page" => "anotherPage.md",
         ],
         format = Documenter.HTML(prettyurls = false)
)
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/forestmod/GenFSM.jl.git",
    devbranch = "main"
)
