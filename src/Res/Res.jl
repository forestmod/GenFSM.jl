module Res

using Proj
import Rasters
import YAML
import DelimitedFiles
import ..ScenarioLoader
SLOAD = ScenarioLoader

include("Pixel.jl")


function load_settings!(settings)
    scenario_path = settings["scenario_path"]

    # Loading res specific settings
    settings["res"] = YAML.load_file(joinpath(scenario_path,"res","settings_res.yaml"))
    # File path in Res settings may contain the {SCENARIO_PATH} placemaker
    SLOAD.recursive_replace!(settings,"\${SCENARIO_PATH}" => scenario_path)

    ft_list = convert(Vector{String},DelimitedFiles.readdlm(settings["res"]["ft_list"],';')[:,1])
    settings["res"]["ft"] = ft_list
    settings["simulation_region"]["nx"] = Int(ceil((settings["simulation_region"]["x_ub"] - settings["simulation_region"]["x_lb"]) / settings["simulation_region"]["xres"])) 
    settings["simulation_region"]["ny"] = Int(ceil((settings["simulation_region"]["y_ub"] - settings["simulation_region"]["y_lb"]) / settings["simulation_region"]["yres"]))
    settings["res"]["nft"] = size(ft_list,1)
    settings["res"]["ndc"] = length(settings["res"]["dc"])

    # Loading res region specific settings
    for reg in settings["res"]["regions"]
        settings["res"][reg] = YAML.load_file(joinpath(scenario_path,"res","settings_res_$(reg).yaml"))
        SLOAD.recursive_replace!(settings["res"][reg],"\${SCENARIO_PATH}" => scenario_path)
    end
end


function init!!(pixels,settings,raster_mask)

    resources_regions = settings["res"]["regions"]

    for (ir, reg) in enumerate(resources_regions)    
        init!!(Val(Symbol(reg)),pixels,settings,raster_mask)
    end

end

end # module Res
