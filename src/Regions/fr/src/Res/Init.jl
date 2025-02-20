

function _init!!(pixels,settings,overal_region_mask)

    println("hello in Res_fr.init!!")
    
    temp_path = joinpath(settings["temp_path"],"res","fr")
    cache_path = joinpath(settings["cache_path"],"res","fr")
    output_path = joinpath(settings["output_path"],"res","fr")
    settings["res"]["fr"]["temp_path"] = temp_path
    settings["res"]["fr"]["cache_path"] = cache_path
    settings["res"]["fr"]["output_path"] = output_path
    isdir(temp_path) || mkpath(temp_path)
    isdir(cache_path) || mkpath(cache_path)
    isdir(output_path) || mkpath(output_path)
    settings["res"]["fr"]["mask"] = get_mask(settings,overal_region_mask)
    mask = Rasters.Raster(settings["res"]["fr"]["mask"])

    get_data!(settings,mask)
    #println(settings)
    
    # Download the data:
    #- DONE administrative for the region
    #- DONE soil 
    #- DONE altimetry DTM
    #- DONE Corine land cover
    #- DONE (no elaboration) to check IGN
    #- TODO Climate, accessibility index
   
end
