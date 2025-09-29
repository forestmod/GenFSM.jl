"""
    module Res

Forest resources module of GenFSM
"""
module Res

using DocStringExtensions
using Proj
import Rasters
import YAML
import DelimitedFiles
import ..ScenarioLoader
SLOAD = ScenarioLoader

include("Pixel.jl")

"""
    load_settings!(settings)

Load the YAML file for the `Res` module and the YAML files for each RES region.

"""
function load_settings!(settings; override=Dict())
    default_scenario_path = settings["default_scenario_path"]
    scenario_path         = settings["scenario_path"]

    # Loading res specific settings
    res_settings = YAML.load_file(joinpath(default_scenario_path,"res","settings_res.yaml"))
    isfile(joinpath(scenario_path,"res","settings_res.yaml")) && SLOAD.override_nested_dict!(res_settings,YAML.load_file(joinpath(scenario_path,"res","settings_res.yaml")))

    settings["res"] = res_settings
    # File path in Res settings may contain the {SCENARIO_PATH} placemaker
    #SLOAD.recursive_replace!(settings,"\${SCENARIO_PATH}" => scenario_path) # no longer needed
    
    # Override Reg settings if called from command line
    # Strict is false because in the override I may have regional stuff that I don't yet have on the main settings that I am overriding
    SLOAD.override_nested_dict!(settings,override, strict=false) 
    #println(scenario_path)
    #println(settings["res"]["ft_list"])
    ft_filepath = isfile(joinpath(scenario_path,settings["res"]["ft_list"])) ? joinpath(scenario_path,settings["res"]["ft_list"]) : joinpath(default_scenario_path,settings["res"]["ft_list"]) 
    ft_list = convert(Vector{String},DelimitedFiles.readdlm(ft_filepath,';')[:,1])
    settings["res"]["ft"] = ft_list
    settings["simulation_region"]["nx"] = Int(ceil((settings["simulation_region"]["x_ub"] - settings["simulation_region"]["x_lb"]) / settings["simulation_region"]["xres"])) 
    settings["simulation_region"]["ny"] = Int(ceil((settings["simulation_region"]["y_ub"] - settings["simulation_region"]["y_lb"]) / settings["simulation_region"]["yres"]))
    settings["res"]["nft"] = size(ft_list,1)
    settings["res"]["ndc"] = length(settings["res"]["dc"])

    # Loading res region specific settings
    for reg in settings["res"]["regions"]
        reg_settings = YAML.load_file(joinpath(default_scenario_path,"res","settings_res_$(reg).yaml"))
        isfile(joinpath(scenario_path,"res","settings_res_$(reg).yaml")) && SLOAD.override_nested_dict!(reg_settings,YAML.load_file(joinpath(scenario_path,"res","settings_res_$(reg).yaml")))
        settings["res"][reg] = reg_settings
        #SLOAD.recursive_replace!(settings["res"][reg],"\${SCENARIO_PATH}" => scenario_path)
    end
end

"""

    $(TYPEDSIGNATURES)

Delegate pixel initialization to each Res region.

Res.init!!() --> init!!(Val(reg_xx)),...) --> Res_xx._init!!()

"""
function init!!(pixels,settings,global_mask)

    resources_regions = settings["res"]["regions"]

    for (ir, reg) in enumerate(resources_regions)    
        init!!(Val(Symbol(reg)),pixels,settings,global_mask)
    end

end

end # module Res
