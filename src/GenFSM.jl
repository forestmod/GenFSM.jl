module GenFSM

#using Pkg 
#Pkg.activate(joinpath(@__DIR__,".."))


import GenFSM_resources as RES
import GenFSM_simrep as REP


export run

include("Utilities.jl")

"""
    run(project,scenario;force_initres_download,force_initres_ml_computation)

Run a scenario of GenFSM

# Parameters
- `project`: project or region name
- `scenario`: name of the scenario within the project
- `force_initres_download`: a dictionary keyed with the resource init regions with the information to what data to force download [def: empty dict]
- `force_initres_ml_train`: a dictionary keyed with the resource init regions with the information to which parts of ML workflow to force retrain [def: empty dict]
- `override`: a vector of pairs (keys=>new value) which eventually overrides the settings as specified in the region/scenario (e.g. for output directories)
"""
function run(project="default",scenario="default";force_initres_download=Dict(),force_initres_ml_train=Dict(),override=[])
    caller_path = pwd()
    settings = REP.load_settings(project,scenario)
    settings["caller_path"] = caller_path
    override_nested_dict!(settings,override) # override settings from command line
    region   = settings["simulation_region"]
    resource_init_regions = settings["resource_init_regions"]
    # first get the raster map and the unistialized pixels
    raster_mask = RES.make_raster_mask(region)
    pixels    = RES.make_pixels(settings["nft"],settings["ndc"],region["nx"],region["ny"])
    #pixels   = RES.init(ft,dc,years,region,init_regions) # return the initialized pixels
    for (ir, reg) in enumerate(resource_init_regions)
        #test = joinpath(@__DIR__,"..","ext","ResourceInit$(reg)Ext.jl")
        #println(test)
        include_expr = Meta.parse("include(joinpath(@__DIR__,\"ext\",\"ResourceInitExt_$(reg).jl\"))") # the parse happens at root level
        eval(include_expr)
        tocall = "resources_init_$(reg)!"

        #toimport = "GenFSM_resource_init_$reg"
        #Pkg.add(toimport)
        #expr = Meta.parse("import $toimport")
        #eval(expr)
        force_download = get(force_initres_download,reg,"")
        force_ml_computation = get(force_initres_ml_computation,reg,"")
        resource_ml_model = Base.invokelatest(getfield(GenFSM, Symbol(tocall)),pixels,raster_mask,settings,force_download,force_ml_computation)
        #RES.set_ml_model(resource_ml_model,ir)
    end
end

function test_path()
    return pwd()
end


end # module GenFSM


