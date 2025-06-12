"""
    GenFSM_simrep

Simulation settings repository for GenFSM models.

This package sole function is [`load_settings`](@ref).

The settings are stored in this repository as yaml files in:

- [project]
  - scenarios
    - [scenario]

Where the default project and the default scenarios ae both named "default".
"""
module ScenarioLoader

using DocStringExtensions
using YAML, DelimitedFiles

import GenFSM

include("Utils.jl")

"""
    load_general_settings(project="default",scenario="default")

Load the settings for a specific project and scenario.
The settings are returned as a dictionary.
"""
function load_general_settings(project="default",scenario="default";override=Dict{Any,Any}())

    # Load default settings
    # Load and override scenario-specific (general, model level) settings
    # Override command line provided settings

    rep_path      = joinpath(@__DIR__,"..","..","repository")
    project_path  = "project_path" in keys(override) ? abspath(override["project_path"]) : joinpath(rep_path,project)
    default_scenario_path = joinpath(project_path,"scenarios","default")
    scenario_path = "scenario_path" in keys(override) ? abspath(override["scenario_path"]) : joinpath(project_path,"scenarios",scenario)
    caller_path   = "caller_path" in keys(override) ? abspath(override["caller_path"]) : pwd()

    # Adding name and path of the region/scenarios to the settings, so that they can be overridden in command line
    settings = Dict{Any,Any}(
        "project"=>project,
        "scenario"=>scenario,
        "default_scenario_path" => default_scenario_path,
        "project_path"=>project_path,
        "scenario_path"=>scenario_path,
        "caller_path"=>caller_path,
    )
    ENV["DEFAULT_SCENARIO_PATH"] = default_scenario_path
    ENV["SCENARIO_PATH"] = scenario_path


    # Loading general scenario settings
    settings      =  merge(settings,YAML.load_file(joinpath(default_scenario_path,"settings.yaml")))
    # Override the defaults with the scenario ones
    override_nested_dict!(settings,YAML.load_file(joinpath(scenario_path,"settings.yaml")))

    # Before the general overriding, checking if the temp_path, cache_path and output_path: (1) has been overriden; (2) they are absolute paths; (3) they exist (or we create them)
    settings["temp_path"] = ("temp_path" in keys(override)) ? ( (isabspath(override["temp_path"])) ? override["temp_path"] : joinpath(caller_path,override["temp_path"])    ) : joinpath(caller_path,project,settings["temp_path"]) 
    settings["cache_path"] = ("cache_path" in keys(override)) ? ( (isabspath(override["cache_path"])) ? override["cache_path"] : joinpath(caller_path,override["cache_path"])    ) : joinpath(caller_path,project,settings["cache_path"]) 
    settings["output_path"] = ("output_path" in keys(override)) ? ( (isabspath(override["output_path"])) ? override["output_path"] : joinpath(caller_path,override["output_path"])    ) : joinpath(caller_path,project,settings["output_path"]) 

    isdir(settings["temp_path"]) || mkpath(settings["temp_path"])
    isdir(settings["cache_path"]) || mkpath(settings["cache_path"])
    isdir(settings["output_path"]) || mkpath(settings["output_path"])
    return settings
end


function load_full_settings(project,scenario;override=Dict())
    # Base.retry_load_extensions() 
    # Load general settings (and early override of paths)
    settings = load_general_settings(project,scenario,override=override)
    # Load settings for the RES module (including regions)
    GenFSM.Res.load_settings!(settings;override=override)

    # Moving verbosity from a string to a enum (integer), so we can say e.g. verbosity <= LOW
    settings["verbosity"] = GenFSM.verbosity_map[settings["verbosity"]]

    # Override all the other settings from command line
    override_nested_dict!(settings,override)
    return settings
end



end # module GenFSM_simrep
