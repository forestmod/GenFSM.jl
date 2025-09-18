"""
    GenFSM module

Main module of the GenFSM package

"""
module GenFSM

using Random
using DocStringExtensions


import DataFrames



export Verbosity, NONE, LOW, STD, HIGH, FULL
export runsim



include("Utils.jl")

global_random_seed::Int = 999 # actual global seed is then read from the YAML file
verbosity::Verbosity = STD

#push!(LOAD_PATH,joinpath(@__DIR__,"ScenarioLoader"))
#import ScenarioLoader as SLOAD
include("ScenarioLoader/ScenarioLoader.jl")
include("Res/Res.jl")

# Regional initializations
include("Regions/fr/src/fr.jl")


#import ScenarioLoader as SLOAD

const SLOAD = ScenarioLoader
const RES   = Res



"""
    run(project,scenario;override)

Run a scenario of GenFSM

# Parameters
- `project`: project or region name
- `scenario`: name of the scenario within the project
- `override`: a dictionary (keys=>new value) which eventually overrides the settings as specified in the region/scenario (e.g. for output directories)
"""
function runsim(project="default",scenario="default";override=Dict{Any,Any}())

    settings = SLOAD.load_full_settings(project,scenario,override=override)
    
    GenFSM.global_random_seed = settings["random_seed"]
    Random.seed!(GenFSM.global_random_seed)
    GenFSM.verbosity = settings["verbosity"]

    region = settings["simulation_region"]
    resources_regions = settings["res"]["regions"]
    # first get the raster map and the unistialized pixels
    # This is the global region mask, and it is always rectangular.
    # Then individual resource or market models can have their own masks
    global_mask = RES.make_global_mask(region)
    pixels      = RES.make_pixels(settings["res"]["nft"],settings["res"]["ndc"],region["nx"],region["ny"])


    RES.init!!(pixels,settings,global_mask)


    # ...
    #rm(settings["temp_path"]; force=true, recursive=true) # TODO: uncomment this
    return settings
end


end # module GenFSM


