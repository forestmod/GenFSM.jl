module GenFSM

#using Pkg 
#Pkg.activate(joinpath(@__DIR__,".."))


import GenFSM_resources as RES
import GenFSM_simrep as REP


export run

include("Utilities.jl")
#include("../ext/ResourceInit.jl")

const DIR_PATH = @__DIR__

"""
    run(project,scenario;override)

Run a scenario of GenFSM

# Parameters
- `project`: project or region name
- `scenario`: name of the scenario within the project
- `override`: a vector of pairs (keys=>new value) which eventually overrides the settings as specified in the region/scenario (e.g. for output directories)
"""
function run(project="default",scenario="default";override=[])
    settings = REP.load_settings(project,scenario,override=override)
    region   = settings["simulation_region"]
    resources_regions = settings["res"]["regions"]
    # first get the raster map and the unistialized pixels
    # This is the global region mask, and it is always rectangular.
    # Then individual resource or market models can have their own masks
    raster_mask = RES.make_raster_mask(region)
    pixels    = RES.make_pixels(settings["res"]["nft"],settings["res"]["ndc"],region["nx"],region["ny"])
    #pixels   = RES.init(ft,dc,years,region,init_regions) # return the initialized pixels

    #println(settings)
    #return nothing

    for (ir, reg) in enumerate(resources_regions)
        
        include_expr = Meta.parse("include(joinpath(DIR_PATH,\"..\",\"ext\",\"ResourceInitExt_$(reg).jl\"))") # the parse happens at root level
        println(include_expr)
        eval(include_expr)
        tocall = "resources_init_$(reg)!"

        #toimport = "GenFSM_resource_init_$reg"
        #Pkg.add(toimport)
        #expr = Meta.parse("import $toimport")
        #eval(expr)
        resource_ml_model = Base.invokelatest(getfield(GenFSM, Symbol(tocall)),pixels,raster_mask,settings)
        #RES.set_ml_model(resource_ml_model,ir)
        #=
        println("debug")
        tocall = "resources_init_$(reg)!"
        resource_ml_model = Base.invokelatest(getfield(GenFSM, Symbol(tocall)),pixels,raster_mask,settings)
        =#
        println("debug ddd") 
        println(resource_ml_model)
    end

    # ...
    rm(settings["temp_path"]; force=true, recursive=true)
    return settings
end


end # module GenFSM


