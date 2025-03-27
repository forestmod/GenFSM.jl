var documenterSearchIndex = {"docs":
[{"location":"anotherPage.html#The-GenFSM-Module","page":"An other page","title":"The GenFSM Module","text":"","category":"section"},{"location":"anotherPage.html","page":"An other page","title":"An other page","text":"GenFSM","category":"page"},{"location":"anotherPage.html#GenFSM","page":"An other page","title":"GenFSM","text":"GenFSM module\n\nMain module of the GenFSM package\n\n\n\n\n\n","category":"module"},{"location":"anotherPage.html#Module-Index","page":"An other page","title":"Module Index","text":"","category":"section"},{"location":"anotherPage.html","page":"An other page","title":"An other page","text":"Modules = [GenFSM]\nOrder   = [:constant, :type, :function, :macro]","category":"page"},{"location":"anotherPage.html#Detailed-API","page":"An other page","title":"Detailed API","text":"","category":"section"},{"location":"anotherPage.html","page":"An other page","title":"An other page","text":"Modules = [GenFSM]\nOrder   = [:constant, :type, :function, :macro]","category":"page"},{"location":"anotherPage.html#GenFSM.runsim","page":"An other page","title":"GenFSM.runsim","text":"run(project,scenario;override)\n\nRun a scenario of GenFSM\n\nParameters\n\nproject: project or region name\nscenario: name of the scenario within the project\noverride: a dictionary (keys=>new value) which eventually overrides the settings as specified in the region/scenario (e.g. for output directories)\n\n\n\n\n\n","category":"function"},{"location":"StyleGuide_templates.html#Style-guide-and-template-for-BetaML-developers","page":"Style guide and template for BetaML developers","title":"Style guide and template for BetaML developers","text":"","category":"section"},{"location":"StyleGuide_templates.html#Master-Style-guide","page":"Style guide and template for BetaML developers","title":"Master Style guide","text":"","category":"section"},{"location":"StyleGuide_templates.html","page":"Style guide and template for BetaML developers","title":"Style guide and template for BetaML developers","text":"The code in GenFSM should follow the official Julia Style Guide.","category":"page"},{"location":"StyleGuide_templates.html#Names-style","page":"Style guide and template for BetaML developers","title":"Names style","text":"","category":"section"},{"location":"StyleGuide_templates.html","page":"Style guide and template for BetaML developers","title":"Style guide and template for BetaML developers","text":"Each file name should start with a capital letter, no spaces allowed (and each file content should start with: \"Part of [GenFSM](https://github.com/forestmod/GenFSM.jl). Licence is MIT.\")\nType names use the so-called \"CamelCase\" convention, where the words are separated by a capital letter rather than _ ,while function names use lower letters only, with words eventually separated (but only when really neeed for readibility) by an _;","category":"page"},{"location":"StyleGuide_templates.html#Docstrings","page":"Style guide and template for BetaML developers","title":"Docstrings","text":"","category":"section"},{"location":"StyleGuide_templates.html","page":"Style guide and template for BetaML developers","title":"Style guide and template for BetaML developers","text":"Please apply the following templates when writing a docstring for GenFSM:","category":"page"},{"location":"StyleGuide_templates.html","page":"Style guide and template for BetaML developers","title":"Style guide and template for BetaML developers","text":"Functions (add @docs if the function is not on the root module level, like for inner constructors, i.e. @docs \"\"\" foo()x ....\"\"\"):","category":"page"},{"location":"StyleGuide_templates.html","page":"Style guide and template for BetaML developers","title":"Style guide and template for BetaML developers","text":"\"\"\"\n    $(TYPEDSIGNATURES)\n\nOne line description\n\n[Further description]\n\n# Parameters:\n\n\n\n# Returns:\n- Elements the funtion need\n\n# Notes:\n- notes\n\n# Example:\n` ` `julia\njulia> [code]\n[output]\n` ` `\n\"\"\"","category":"page"},{"location":"StyleGuide_templates.html","page":"Style guide and template for BetaML developers","title":"Style guide and template for BetaML developers","text":"Structs","category":"page"},{"location":"StyleGuide_templates.html","page":"Style guide and template for BetaML developers","title":"Style guide and template for BetaML developers","text":"\"\"\"\n$(TYPEDEF)\n\nOne line description\n\n[Further description]\n\n# Fields: (if relevant)\n$(TYPEDFIELDS)\n\n# Notes:\n\n# Example:\n` ` `julia\njulia> [code]\n[output]\n` ` `\n\n\"\"\"","category":"page"},{"location":"StyleGuide_templates.html","page":"Style guide and template for BetaML developers","title":"Style guide and template for BetaML developers","text":"Enums:","category":"page"},{"location":"StyleGuide_templates.html","page":"Style guide and template for BetaML developers","title":"Style guide and template for BetaML developers","text":"\"\"\"\n$(TYPEDEF)\n\nOne line description\n\n[Further description]\n\n\n# Notes:\n\n\"\"\"","category":"page"},{"location":"StyleGuide_templates.html","page":"Style guide and template for BetaML developers","title":"Style guide and template for BetaML developers","text":"Constants","category":"page"},{"location":"StyleGuide_templates.html","page":"Style guide and template for BetaML developers","title":"Style guide and template for BetaML developers","text":"\"\"\"\n[4 spaces] [Constant name]\n\nOne line description\n\n[Further description]\n\n\n# Notes:\n\n\"\"\"","category":"page"},{"location":"StyleGuide_templates.html","page":"Style guide and template for BetaML developers","title":"Style guide and template for BetaML developers","text":"Modules","category":"page"},{"location":"StyleGuide_templates.html","page":"Style guide and template for BetaML developers","title":"Style guide and template for BetaML developers","text":"\"\"\"\n[4 spaces] [Module name]\n\nOne line description\n\nDetailed description on the module objectives, content and organisation\n\n\"\"\"","category":"page"},{"location":"StyleGuide_templates.html#Internal-links","page":"Style guide and template for BetaML developers","title":"Internal links","text":"","category":"section"},{"location":"StyleGuide_templates.html","page":"Style guide and template for BetaML developers","title":"Style guide and template for BetaML developers","text":"To refer to a documented object: [`NAME`](@ref) or [`NAME`](@ref manual_id). In particular for internal links use [`?NAME`](@ref ?NAME)","category":"page"},{"location":"StyleGuide_templates.html","page":"Style guide and template for BetaML developers","title":"Style guide and template for BetaML developers","text":"To create an id manually: [Title](@id manual_id)","category":"page"},{"location":"StyleGuide_templates.html#Data-organisation","page":"Style guide and template for BetaML developers","title":"Data organisation","text":"","category":"section"},{"location":"StyleGuide_templates.html","page":"Style guide and template for BetaML developers","title":"Style guide and template for BetaML developers","text":"While some functions provide a dims parameter, most GenFSM algorithms expect the input data layout with observations organised by rows and fields/features by columns.\nWhile some algorithms accept as input DataFrames, the usage of standard arrays is encourages (if the data is passed to the function as dataframe, it may be converted to standard arrays somewhere inside inner loops, leading to great inefficiencies).","category":"page"},{"location":"Structure_of_model.html#Model-code-structure","page":"Model code structure","title":"Model code structure","text":"","category":"section"},{"location":"Structure_of_model.html","page":"Model code structure","title":"Model code structure","text":"The following is a simplification of how the model interact between the various module, including regional specific ones:","category":"page"},{"location":"Structure_of_model.html","page":"Model code structure","title":"Model code structure","text":"module GenFSM\n\nexport run_model\n\nmodule Res\n\nfunction init!!()\n    res_regions = [\"fr\",\"uk\"]\n    println(\"Res init..\")\n    for reg in res_regions   \n        init!!(Val(Symbol(reg)),1)\n    end\nend\n\nend # end module Res\n\nmodule Mkt\n\nfunction init!!()\n    mkt_regions = [\"fr\",\"de\"]\n    println(\"Mkt init..\")\n    for reg in mkt_regions   \n        init!!(Val(Symbol(reg)),1)\n    end\nend\n\nend # end module mkt\n\nmodule Res_fr\nimport ..GenFSM.Res\n\nfunction _init!!(x)\n    println(\"Res fr\")\nend\nfunction Res.init!!(::Val{:fr},x)\n    Res_fr._init!!(x)\nend\n\nend\n\nmodule Res_uk\nimport ..GenFSM.Res\n\nfunction _init!!(x)\n    println(\"Res uk\")\nend\nfunction Res.init!!(::Val{:uk},x)\n    Res_uk._init!!(x)\nend\n\nend\n\nmodule Mkt_fr\nimport ..GenFSM.Mkt\n\nfunction _init!!(x)\n    println(\"mkt fr\")\nend\nfunction Mkt.init!!(::Val{:fr},x)\n    Mkt_fr._init!!(x)\nend\n\nend\n\nmodule Mkt_de\nimport ..GenFSM.Mkt\n\nfunction _init!!(x)\n    println(\"mkt de\")\nend\nfunction Mkt.init!!(::Val{:de},x)\n    Mkt_de._init!!(x)\nend\n\nend # end Mkt_de\n\nfunction run_model()\n    Res.init!!()\n    Mkt.init!!()\nend\n\nend # end master module GenFSM\n\n\nGenFSM.run_model()","category":"page"},{"location":"index.html#GenFSM.jl","page":"Index","title":"GenFSM.jl","text":"","category":"section"},{"location":"index.html","page":"Index","title":"Index","text":"Documentation for GenFSM.jl","category":"page"}]
}
