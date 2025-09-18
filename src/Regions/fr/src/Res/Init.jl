
"""
    $(TYPEDSIGNATURES)

Initialization function of the `Res` module for the French region


    
"""
function _init!!(pixels,settings,global_mask)

    println("Hello in Res_fr.init!!")
    
    temp_path = joinpath(settings["temp_path"],"res","fr")
    cache_path = joinpath(settings["cache_path"],"res","fr")
    output_path = joinpath(settings["output_path"],"res","fr")
    settings["res"]["fr"]["temp_path"] = temp_path
    settings["res"]["fr"]["cache_path"] = cache_path
    settings["res"]["fr"]["output_path"] = output_path
    isdir(temp_path) || mkpath(temp_path)
    isdir(cache_path) || mkpath(cache_path)
    isdir(output_path) || mkpath(output_path)
    settings["res"]["fr"]["mask"] = make_reg_res_mask(settings,global_mask)
    reg_res_mask = Rasters.Raster(settings["res"]["fr"]["mask"]) # still the inner matrix is an Union{Missing,Int62}
    reg_res_mask = map(i -> i == 0 ? 0 : 1, reg_res_mask[:,:])
    get_data!(settings,reg_res_mask)
    #println(settings)
    
    # Download the data:
    #- DONE administrative for the region
    #- DONE soil 
    #- DONE altimetry DTM
    #- DONE Corine land cover
    #- DONE (no elaboration) to check IGN
    #- DONE Climate
    # Maybe TODO: accessibility index
   
    prepare_data!(settings,reg_res_mask)
    
end
