# See also Copernicus or CHELSA for climate:
# https://cds.climate.copernicus.eu/cdsapp#!/search?type=dataset
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=overview
# https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/EUR11/documentation/CHELSA_EUR11_technical_documentation.pdf


function get_data!(settings,mask)
    @info  "Getting data for the resource initialization for France..."
    input_rasters = Dict{Any,Any}()
    @info  "- getting dtm..."
    input_rasters["dtm"] = get_dtm(settings,mask)
    #input_rasters = merge(input_rasters,get_dtm(settings,mask))
    @info  "- getting soil data..."
    input_rasters["soil"] = get_soil_data(settings,mask)
    #input_rasters = merge(input_rasters,get_soil_data(settings,mask))
    @info  "- getting clc..."
    input_rasters["clc"] = get_clc(settings,mask)
    #input_rasters = merge(input_rasters,get_clc(settings,mask))
    @info  "- getting climate (historic or future depending on scenario)..."
    input_rasters["clim"] = get_climate_data!(settings,mask)
    settings["res"]["fr"]["input_rasters"] = input_rasters
    @info  "- getting CO2 atmospheric concentration (historic or future depending on scenario)..."
    co2_conc_h, co2_conc_f = get_co2_concentration(settings,mask)
    settings["res"]["fr"]["co2_conc_h"] = co2_conc_h; settings["res"]["fr"]["co2_conc_f"] = co2_conc_f

    @info  "- getting inventory data..."
    inv_data = get_inventory_data(settings,mask)
    settings["res"]["fr"]["inv_data"] = inv_data

    @info  "- DONE getting data"

end

# ------------------------------------------------------------------------------
# Specific download functions

function get_mask(settings,mask)
    force          = "adm_borders" in settings["res"]["fr"]["force_download"]
    verbosity      = settings["verbosity"]
    verbose        = verbosity >= HIGH
    admin_borders_input_crs = settings["res"]["fr"]["data_sources"]["admin_borders_input_crs"]

    adm_borders_path      = joinpath(settings["res"]["fr"]["cache_path"],"adm_borders")
    isdir(adm_borders_path) || mkpath(adm_borders_path)
    adm_borders_dl_path   = joinpath(adm_borders_path,"downloaded")
    isdir(adm_borders_dl_path) || mkpath(adm_borders_dl_path)

    mask_destpath         = joinpath(adm_borders_path,"mask.tif")

    (isfile(mask_destpath) && (!force)) && return mask_destpath

    for f in settings["res"]["fr"]["data_sources"]["admin_borders_sources"]
        fname     = basename(split(f,"?")[1])
        dest_name = joinpath(adm_borders_dl_path,fname)
        Downloads.download(f,dest_name, verbose=verbose)
    end
    dir_files = readdir(adm_borders_dl_path)
    shp_file = dir_files[findfirst(x->match(r".*\.shp$", x) !== nothing, dir_files)]

    reg_borders = Shapefile.Handle(joinpath(adm_borders_dl_path,shp_file)).shapes

    reg_raster = Rasters.rasterize(last, reg_borders; res=0.01, missingval=0, fill=1, progress=true)
    shp_crs           = convert(Rasters.WellKnownText, Rasters.EPSG(admin_borders_input_crs))
    reg_raster = Rasters.setcrs(reg_raster, shp_crs)
    reg_raster  = Rasters.reverse(reg_raster;dims=Rasters.Y )
    #Rasters.metadata(reg_raster)["missed_geometries"]

    resampled_raster = Rasters.resample(reg_raster,to=mask,method=:average)
    write(mask_destpath, resampled_raster ,force=true)
    rm(adm_borders_dl_path,recursive=true)
    return mask_destpath
end


"""
   get_dtm(settings,mask)

Download and resample to mask dtm and related variables (slope, aspect, 
Terrain Ruggedness Index)

"""
function get_dtm(settings,mask)
    data_path        = joinpath(settings["res"]["fr"]["cache_path"],"dtm")
    isdir(data_path) || mkpath(data_path)
    force            = "dtm" in settings["res"]["fr"]["force_download"]
    url              = settings["res"]["fr"]["data_sources"]["dtm_url"]
    to               = mask
    filename         = basename(url)
    verbosity        = settings["verbosity"]
    verbose          = verbosity in [HIGH, FULL] 
    zip_destpath     = joinpath(data_path,filename)
    tif_destpath     = replace(zip_destpath,".zip" => ".tif")
    tif_destpath_reg = replace(tif_destpath,".tif" => "_reg.tif")
    slope_destpath   = joinpath(data_path,"slope.tif")
    aspect_destpath  = joinpath(data_path,"aspect.tif")
    tri_destpath     = joinpath(data_path,"tri.tif")
    dtm_var          = DataStructures.OrderedDict(
        "dtm"=>tif_destpath_reg,
        "slope"=>slope_destpath,
        "aspect"=>aspect_destpath,
        "tri"=>tri_destpath,
        )
    if ( isfile(tif_destpath_reg) && (!force))
        return dtm_var
    end
    verbose && @info "Downloading dtm file..."
    Downloads.download(url,zip_destpath, verbose=verbose)
    unzip(zip_destpath,data_path)
    sr            = settings["simulation_region"]
    crs           = convert(Rasters.WellKnownText, Rasters.EPSG(sr["cres_epsg_id"]))
    lon, lat      = Rasters.X(sr["x_lb"]:(sr["xres"]/4):sr["x_ub"]), Rasters.Y(sr["y_lb"]:(sr["yres"]/4):sr["y_ub"])
    mask_hr       = Rasters.Raster(zeros(Float32,lon,lat),crs=crs)
    mask_hr       = Rasters.reverse(mask_hr;dims=Rasters.Y )
    dtm_w         = Rasters.Raster(tif_destpath)
    dtm_reg_hr    = Rasters.resample(dtm_w,to=mask_hr,method=:average)
    dtm_reg       = Rasters.resample(dtm_w,to=to,method=:average)
    M             = convert(Matrix{Float32},dtm_reg_hr[:,:])
    slope         = Geomorphometry.slope(M)
    aspect        = Geomorphometry.aspect(M)
    tri           = Geomorphometry.TRI(M)
    slope_r       = similar(mask_hr)
    aspect_r      = similar(mask_hr)
    tri_r         = similar(mask_hr)
    slope_r[:,:]  = slope
    aspect_r[:,:] = aspect
    tri_r[:,:]    = tri
    slope_r_lr    = Rasters.resample(slope_r,to=to,method=:average)
    aspect_r_lr   = Rasters.resample(aspect_r,to=to,method=:average)
    tri_r_lr      = Rasters.resample(tri_r,to=to,method=:average)
    write(tif_destpath_reg, dtm_reg,force=true)
    write(slope_destpath, slope_r_lr,force=true)
    write(aspect_destpath, aspect_r_lr,force=true)
    write(tri_destpath, tri_r_lr,force=true)
    rm(zip_destpath)
    rm(tif_destpath)
    return dtm_var
end

function get_soil_data(settings,mask)
    #=
    All soil datasets:
    https://esdac.jrc.ec.europa.eu/resource-type/european-soil-database-soil-properties

    Soil_physical 
    https://esdac.jrc.ec.europa.eu/content/topsoil-physical-properties-europe-based-lucas-topsoil-data
    https://www.sciencedirect.com/science/article/pii/S0016706115300173

    Soil chemistry
    https://esdac.jrc.ec.europa.eu/content/chemical-properties-european-scale-based-lucas-topsoil-data
    https://www.sciencedirect.com/science/article/pii/S0016706119304768

    Soil other derived data
    https://esdac.jrc.ec.europa.eu/content/european-soil-database-derived-data
    Hiederer, R. 2013. Mapping Soil Properties for Europe - Spatial Representation of Soil Database Attributes. Luxembourg: Publications Office of the European Union - 2013 - 47pp. EUR26082EN Scientific and Technical Research series, ISSN 1831-9424, doi:10.2788/94128
    Hiederer, R. 2013. Mapping Soil Typologies - Spatial Decision Support Applied to European Soil Database. Luxembourg: Publications Office of the European Union - 2013 - 147pp. EUR25932EN Scientific and Technical Research series, ISSN 1831-9424, doi:10.2788/8728

    Look also this: https://www.isric.org/explore/soilgrids
    And lucas: calcares

    =#

    soil_path      = joinpath(settings["res"]["fr"]["cache_path"],"soil")
    isdir(soil_path) || mkpath(soil_path)
    soil_ph_url    = settings["res"]["fr"]["data_sources"]["soil_ph_url"]
    soil_chem_url  = settings["res"]["fr"]["data_sources"]["soil_chem_url"]
    soil_oth_url   = settings["res"]["fr"]["data_sources"]["soil_oth_url"]
    soil_ph_vars   = settings["res"]["fr"]["data_sources"]["soil_ph_vars"]
    soil_chem_vars = settings["res"]["fr"]["data_sources"]["soil_chem_vars"]
    soil_oth_vars  = settings["res"]["fr"]["data_sources"]["soil_oth_vars"]
    force          = "soil" in settings["res"]["fr"]["force_download"]
    verbosity      = settings["verbosity"]
    verbose        = verbosity in [HIGH, FULL] 
    to             = mask
    soil_vars      = DataStructures.OrderedDict{String,String}()
    n_soil_ph_vars = length(soil_ph_vars)
    soil_texture_n_classes = settings["res"]["fr"]["data_sources"]["soil_texture_n_classes"]
    crs = convert(Rasters.WellKnownText, Rasters.EPSG(settings["simulation_region"]["cres_epsg_id"]))

    # Soil and chemistry variables....
    for (i,var) in enumerate(vcat(soil_ph_vars,soil_chem_vars))
        if i <= n_soil_ph_vars
            urlname   = replace(soil_ph_url,"\${VAR}" => var)
        else
            urlname   = replace(soil_chem_url,"\${VAR}" => var)
        end
        zipname   = joinpath(soil_path,"$(var).zip")
        zipfolder = joinpath(soil_path,var)
        final_file = joinpath(soil_path,"$(var)_reg.tif") 
        if (ispath(final_file) && !force )
            soil_vars[var] = final_file
            continue
        end
        Downloads.download(urlname,zipname, verbose=verbose)
        unzip(zipname,zipfolder)
        dir_files = readdir(zipfolder)
        tif_file = dir_files[findfirst(x->match(r".*\.tif$", x) !== nothing, dir_files)]
        saved_file = joinpath(zipfolder,tif_file)
        orig_raster = Rasters.Raster(saved_file)
        resampled_raster = Rasters.resample(orig_raster,to=to,method=:average)
        write(final_file, resampled_raster, force=true)
        soil_vars[var] = final_file
        rm(zipname)
        rm(zipfolder,recursive=true)
    end

    # Special: for TextureUSDA we extract the class as boolean values, as this is categorical
    soil_texture_classes = 1:soil_texture_n_classes
    texture_filename = joinpath(soil_path,"TextureUSDA_reg.tif")
    texture_classes = expand_classes(texture_filename,soil_texture_classes;verbose=verbose,force=force,to=nothing)
    delete!(soil_vars,"TextureUSDA")
    soil_vars = DataStructures.OrderedDict(soil_vars..., texture_classes...)

    # Other soil variables
    urlname   = soil_oth_url
    zipname   = joinpath(soil_path,basename(urlname))
    zipfolder = joinpath(soil_path, split(basename(urlname),".")[1])
    finalfolder = joinpath(soil_path,"soil_oth_vars")
    if (ispath(finalfolder) && (!force) )
        return soil_vars
    end
    ispath(finalfolder) || mkpath(finalfolder)
    Downloads.download(urlname,zipname, verbose=verbose)
    unzip(zipname,zipfolder)
    for var in soil_oth_vars
        saved_file = joinpath(zipfolder,"$(var).rst")
        final_file = joinpath(finalfolder,"$(var)_reg.tif")
        orig_raster = Rasters.Raster(saved_file,crs=crs) # ,missingval=0.0f0)  missing value is zero, but zero is used also as true value! 
        resampled_raster = Rasters.resample(orig_raster,to=to,method=:average)
        write(final_file, resampled_raster,force=true)
        soil_vars[var] = final_file
    end
    rm(zipname)
    rm(zipfolder,recursive=true)
    return soil_vars
end

function get_clc(settings,mask)
    clc_dirpath   = joinpath(settings["res"]["fr"]["cache_path"],"clc")
    clc_dldirpath = joinpath(settings["res"]["fr"]["temp_path"],"clc")
    clc_dlpath    = joinpath(clc_dldirpath,"clc.gpkg")
    clc_bfor_path = joinpath(clc_dirpath,"clc_bfor.tif")
    clc_cfor_path = joinpath(clc_dirpath,"clc_cfor.tif")
    clc_mfor_path = joinpath(clc_dirpath,"clc_mfor.tif")
    clc_vars = Dict(
        "clc_bfor" => clc_bfor_path,
        "clc_cfor" => clc_cfor_path,
        "clc_mfor" => clc_mfor_path,
    )
    clc_url = settings["res"]["fr"]["data_sources"]["clc_url"]
    force          = "clc" in settings["res"]["fr"]["force_download"]
    (isdir(clc_dirpath) && (!force) ) && return clc_vars
    isdir(clc_dirpath) || mkpath(clc_dirpath)
    isdir(clc_dldirpath) || mkpath(clc_dldirpath)
    Downloads.download(clc_url,clc_dlpath)
    clc_df  = GeoDataFrames.read(clc_dlpath)
    DataFrames.rename!(clc_df,"Shape" => "geometry")
    clc_bfor = clc_df[clc_df.Code_18 .== "311",:]
    clc_cfor = clc_df[clc_df.Code_18 .== "312",:]
    clc_mfor = clc_df[clc_df.Code_18 .== "313",:]
    Logging.with_logger(Logging.NullLogger()) do
        clc_bfor_share = Rasters.coverage(clc_bfor; to=mask)
        clc_cfor_share = Rasters.coverage(clc_cfor; to=mask)
        clc_mfor_share = Rasters.coverage(clc_mfor; to=mask)
        write(clc_bfor_path, clc_bfor_share ,force=true)
        write(clc_cfor_path, clc_cfor_share ,force=true)
        write(clc_mfor_path, clc_mfor_share ,force=true)
    end
    rm(clc_dldirpath, recursive=true)
    return clc_vars
end

"""
   get_climate_data!(settings,mask)

Download the historical and (eventually) the future climate data

Look also:
# AGera5 0.1 degrees
# Agri4Cas: 0.3 degrees
"""
function get_climate_data!(settings,mask)

    verbosity      = settings["verbosity"]
    verbose        = verbosity in [HIGH, FULL] 
    
    # Remember that scenario is already embedded in the temp path, but the cache path shares several scenarios
    cl_h_temppath        = joinpath(settings["res"]["fr"]["temp_path"],"clim_historical")
    cl_h_cachepath       = joinpath(settings["res"]["fr"]["cache_path"],"clim_historical")
    cl_f_temppath        = joinpath(settings["res"]["fr"]["temp_path"],"clim_future")
    cl_f_cachepath       = joinpath(settings["res"]["fr"]["cache_path"],"clim_future",settings["scenario"])

    clim_settings = settings["res"]["fr"]["data_sources"]["clim"]
    force_h  = "clim_h" in settings["res"]["fr"]["force_download"]
    force_f  = "clim_f" in settings["res"]["fr"]["force_download"]  
    freeze =  clim_settings["fixed_climate"]
    vars   =  clim_settings["vars"]
    hyears =  clim_settings["hist_years"]
    fyears =  parse(Int64,clim_settings["fefps"][1:4]):parse(Int64,clim_settings["fefpe"][1:4])
    fyears =  maximum(hyears)+1:fyears[end]
    settings["res"]["fr"]["data_sources"]["clim"]["fut_years"] = fyears # updating the settings with the future years
    mmonths = [m<10 ? "0$m" : "$m" for m in 1:12]
    tr_fs_h = clim_settings["transformations_h"]
    tr_fs_f = clim_settings["transformations_f"]

    hclim_data      = DataStructures.OrderedDict{Tuple,String}()
    fclim_data      = DataStructures.OrderedDict{Tuple,String}()

    [hclim_data[(var,m,y)] = joinpath(cl_h_cachepath,"var_$(var)_$(m)_$(y).tif") for var in vars, m in 1:12, y in hyears]

    [fclim_data[(var,m,y)] = joinpath(cl_f_cachepath,"var_$(var)_$(m)_$(y).tif") for var in vars, m in 1:12, y in fyears]

    # Download the historical dataset...
    if !isdir(cl_h_cachepath) || force_h == true

        isdir(cl_h_temppath)  || mkpath(cl_h_temppath)
        isdir(cl_h_cachepath) || mkpath(cl_h_cachepath)
        @info "Downloading raw historical climatic data..."

        hist_base_url = clim_settings["hist_base_url"]

        
        for v in vars, m in 1:12
            tr_f_str = (tr_fs_h == nothing) ? "identity" : get(tr_fs_h,v,"identity")
            tr_f = eval(Meta.parse(tr_f_str))
            for y in hyears
                mm  = mmonths[m]
                url = replace(hist_base_url, "\${VAR}" => v, "\${MMONTH}" => mm, "\${YEAR}" => y )
                temp_path  = joinpath(cl_h_temppath,"var_$(v)_$(m)_$(y).tif")
                final_path = joinpath(cl_h_cachepath,"var_$(v)_$(m)_$(y).tif")
                Downloads.download(url,temp_path, verbose=verbose)
                orig_raster    = Rasters.Raster(temp_path)
                orig_raster    = invokelatest(Rasters.modify,tr_f,orig_raster) # transformation to align measure unit or to correct bias
                resampled_raster = Rasters.resample(orig_raster,to=mask,method=:average)
                write(final_path, resampled_raster, force=true)
                rm(temp_path) # I should not need to use neither force nor recursive
            end
        end
        
    end

    # Downloading scenario-based future climatic data
    (force_f == false && isdir(cl_f_cachepath)) &&  return Dict("historical"=>hclim_data, "future"=>fclim_data)

    isdir(cl_f_temppath)  || mkpath(cl_f_temppath)
    isdir(cl_f_cachepath) || mkpath(cl_f_cachepath)

    if clim_settings["fixed_climate"]
        @info "Copying historical clim data as future data..."
        # Setting as future weather the last 10 years of obs weather, cyclically
        # The cycle reverses direction at each last_obs_years step
        last_obs_years = hyears[max(1,length(hyears)-10):end]
        deltay = -1
        ihy    = length(last_obs_years)
        for (iy, y) in enumerate(fyears)
            hy = last_obs_years[ihy]
            for v in vars, m in 1:12
                source_path = joinpath(cl_h_cachepath,"var_$(v)_$(m)_$(hy).tif")
                source_rel_path = joinpath("..","..","clim_historical","var_$(v)_$(m)_$(hy).tif")
                dest_path   = joinpath(cl_f_cachepath,"var_$(v)_$(m)_$(y).tif")
                isfile(dest_path) && rm(dest_path)
                try
                    if Sys.iswindows() && Sys.windows_version() < Sys.WINDOWS_VISTA_VER
                        cp(source_path,dest_path)
                    else
                        symlink(source_rel_path, dest_path)
                    end
                catch
                    cp(source_path,dest_path)
                end
            end
            # updatting index for next yearly loop
            (ihy == 1) && (deltay = +1)
            (ihy == length(last_obs_years)) && (deltay = -1)
            ihy = ihy + deltay
        end
        return Dict("historical"=>hclim_data, "future"=>fclim_data)
    end

    # Not in a freezing scenario, let's download scenario-specific future climatic data...
    @info "Downloading climatic future data..."

    fut_base_path = clim_settings["fut_base_path"]
    projid    = settings["simulation_region"]["cres_epsg_id"]
    x_lb_proj = settings["simulation_region"]["x_lb"]
    y_lb_proj = settings["simulation_region"]["y_lb"]
    x_ub_proj = settings["simulation_region"]["x_ub"]
    y_ub_proj = settings["simulation_region"]["y_ub"]
    inv_trans = Proj.Transformation("EPSG:$(projid)","EPSG:4326",always_xy=true)

    # Trasformation in geographical coordinates for the area to download
    x_lb1,_ = inv_trans(x_lb_proj,y_lb_proj)
    x_lb2,_ = inv_trans(x_lb_proj,y_ub_proj)
    x_lb    = min(x_lb1,x_lb2) - 0.05
    _,y_lb1 = inv_trans(x_lb_proj,y_lb_proj)
    _,y_lb2 = inv_trans(x_ub_proj,y_lb_proj)
    y_lb    = min(y_lb1,y_lb2) - 0.05
    x_ub1,_ = inv_trans(x_ub_proj,y_lb_proj)
    x_ub2,_ = inv_trans(x_ub_proj,y_ub_proj)
    x_ub    = max(x_ub1,x_ub2) + 0.05
    _,y_ub1 = inv_trans(x_lb_proj,y_ub_proj)
    _,y_ub2 = inv_trans(x_ub_proj,y_ub_proj)
    y_ub    = max(y_ub1,y_ub2) + 0.05

    for fy in fyears
        ## Download year data
        cl_f_temppath_y = joinpath(cl_f_temppath,"$(fy)/",)
        isdir(cl_f_temppath_y)  || mkpath(cl_f_temppath_y)
    
        CondaPkg.withenv() do
            # All this because it doesn't work with just the function call
            python = CondaPkg.which("python")
            chelsa_script = joinpath(cl_f_temppath,"00_chelsa_dl_cmip6_$(fy).py")
    
            write(chelsa_script,
"""
from chelsa_cmip6.GetClim import chelsa_cmip6

chelsa_cmip6(
    activity_id='$(clim_settings["activity_id"])', 
    table_id='$(clim_settings["table_id"])', 
    experiment_id='$(clim_settings["experiment_id"])', 
    institution_id='$(clim_settings["institution_id"])', 
    source_id='$(clim_settings["source_id"])', 
    member_id='$(clim_settings["member_id"])', 
    refps='$(clim_settings["refps"])', 
    refpe='$(clim_settings["refpe"])', 
    fefps='$(fy)$(clim_settings["fefps"][5:end])', 
    fefpe='$(fy)$(clim_settings["fefpe"][5:end])', 
    xmin=$(x_lb), 
    xmax=$(x_ub),
    ymin=$(y_lb), 
    ymax=$(y_ub),
    output='$(cl_f_temppath_y)',
    use_esgf=False,
    bio=False
)
""")
           run(`$python $(chelsa_script)`)
        end

        for v in vars
            tr_f_str = (tr_fs_f == nothing) ? "identity" : get(tr_fs_f,v,"identity")
            tr_f = eval(Meta.parse(tr_f_str))
            saved_path =  replace(fut_base_path,
                "\${INSTITUTION_ID}" => clim_settings["institution_id"],
                "\${SOURCE_ID}" => clim_settings["source_id"],
                "\${VAR}" => v,
                "\${EXPERIMENT_ID}" => clim_settings["experiment_id"],
                "\${MEMBER_ID}" => clim_settings["member_id"],
                "\${FEFPS}" => "$(fy)$(clim_settings["fefps"][5:end])",
                "\${FEFPE}" => "$(fy)$(clim_settings["fefpe"][5:end])",
                )
            saved_fullpath = joinpath(cl_f_temppath_y,saved_path)
            var_monthly  = Rasters.Raster(saved_fullpath; name=v)
            var_monthly = invokelatest(Rasters.modify,tr_f,var_monthly)  # transformation for bias or align measure unit
            var_monthly = Rasters.reverse(var_monthly, dims=Rasters.Y) # the nc data is not in the North-first format
            resampled_raster = Rasters.resample(var_monthly,to=mask,method=:average)
            for m in 1:12 
                final_path = joinpath(cl_f_cachepath,"var_$(v)_$(m)_$(fy).tif")
                write(final_path, resampled_raster[:,:,m], force=true)
            end 
        end
        rm(cl_f_temppath_y, recursive=true)
    end # end of each fyear
    
    return Dict("historical"=>hclim_data, "future"=>fclim_data)


end

"""
   get_co2_concentration(settings,mask)

Download and retrieve the yearly average CO2 concentration in the atmosphere, historical and scenario-dependant future one.
"""
function get_co2_concentration(settings,mask)

    verbosity      = settings["verbosity"]
    verbose        = verbosity in [GenFSM.HIGH, GenFSM.FULL] 

    # Remember that scenario is already embedded in the temp path, but the cache path shares several scenarios
    co2_h_destfolder = joinpath(settings["res"]["fr"]["cache_path"],"co2_historical")
    co2_f_destfolder = joinpath(settings["res"]["fr"]["cache_path"],"co2_future",settings["scenario"])
    isdir(co2_h_destfolder) || mkpath(co2_h_destfolder)
    isdir(co2_f_destfolder) || mkpath(co2_f_destfolder)
    co2_h_destpath       = joinpath(co2_h_destfolder,"co2_conc.csv")
    co2_f_destpath       = joinpath(co2_f_destfolder,"co2_conc.csv")
    clim_settings = settings["res"]["fr"]["data_sources"]["clim"]
    force_h  = "co2_conc_h" in settings["res"]["fr"]["force_download"]
    force_f  = "co2_conc_f" in settings["res"]["fr"]["force_download"]  
    freeze =  clim_settings["fixed_climate"] # boolean
    hyears =  clim_settings["hist_years"]
    fyears =  clim_settings["fut_years"] # added in get_climate_data!

    if !force_h && !force_f && ispath(co2_h_destpath) && ispath(co2_f_destpath)
        return (co2_h_destpath, co2_f_destpath)
    end

    # ok, we need to do somehting, let's download the file and get the data
    co2_conc_url = settings["res"]["fr"]["data_sources"]["co2_conc_url"]
    co2_conc_url = replace(co2_conc_url,"\${SSP_SCENARIO}" => uppercase(clim_settings["experiment_id"]))  # eg. "ssp585" or "ssp126

    data = @pipe HTTP.get(co2_conc_url).body |>
        CSV.File(_,delim=' ',header=false, ignorerepeated=true, skipto=3) |> DataFrames.DataFrame
    co2_h = data[in.(data[:,1],Ref(hyears)),2]
    if freeze
        co2_f = fill(co2_h[end], length(fyears))
    else
        co2_f = data[in.(data[:,1],Ref(fyears)),2]
    end

    if force_h || ! ispath(co2_h_destpath) 
    CSV.write(co2_h_destpath, DataFrames.DataFrame(years=hyears,co2_conc=co2_h))
    end
    if force_f || ! ispath(co2_f_destpath) 
    CSV.write(co2_f_destpath, DataFrames.DataFrame(years=fyears,co2_conc=co2_f))
    end
    return (co2_h_destpath, co2_f_destpath)
end


"""
   get_inventory_data(settings,mask)

Download the forest inventory data and change the (X,Y) coordinated of the points to the CRS used in the model.
Returns a dicionary with the file paths (no further elaborated, still in csv format)   
"""
function get_inventory_data(settings,mask)
    forinv_url       = settings["res"]["fr"]["data_sources"]["forest_inventory_url"]
    force            = "forinv" in settings["res"]["fr"]["force_download"]
    forest_inventory_crs = settings["res"]["fr"]["data_sources"]["forest_inventory_cres_epsg_id"]
    forinv_dirpath   = joinpath(settings["res"]["fr"]["cache_path"],"forinv")
    forinv_dldirpath = joinpath(settings["res"]["fr"]["temp_path"],"forinv")
    forinv_dlpath    = joinpath(forinv_dldirpath,basename(forinv_url))
    forinv_unzippeddir = joinpath(forinv_dldirpath,"data")
    forinv_data = Dict(
        "points"          =>  joinpath(forinv_dirpath,"PLACETTE.csv"),
        "trees"           =>  joinpath(forinv_dirpath,"ARBRE.csv"),
        "death_trees"     =>  joinpath(forinv_dirpath,"BOIS_MORT.csv"),
        "tree_cover"      =>  joinpath(forinv_dirpath,"COUVERT.csv"),
        "points_toposoil" =>  joinpath(forinv_dirpath,"ECOLOGIE.csv"),
        "species"         =>  joinpath(forinv_dirpath,"FLORE.csv"),
        "habitat"         =>  joinpath(forinv_dirpath,"HABITAT.csv"),
    )
    forinv_meta = Dict(
        "species_latin_name"        => joinpath(forinv_dirpath,"espar-cdref13.csv"),
        "vars_availability_by_year" => joinpath(forinv_dirpath,"summary_data.csv"), # recapitulatif_donnees.csv is the French version
        "metadata"                  => joinpath(forinv_dirpath,"metadata.csv"), # metadonnees.csv is the French version    
    )
    forinv_vars = Dict("data" => forinv_data, "meta"=>forinv_meta)
    (isdir(forinv_dirpath) && (!force) ) && return forinv_vars
    isdir(forinv_dirpath) || mkpath(forinv_dirpath)
    isdir(forinv_dldirpath) || mkpath(forinv_dldirpath)
    Downloads.download(forinv_url,forinv_dlpath)
    unzip(forinv_dlpath,forinv_unzippeddir)
    mv(forinv_unzippeddir,forinv_dirpath,force=true)
    points = CSV.read(forinv_data["points"],DataFrames.DataFrame)
    trans = Proj.Transformation("EPSG:$(forest_inventory_crs)", "EPSG:$(settings["simulation_region"]["cres_epsg_id"])", always_xy=true)   
    for p in eachrow(points)
      (X,Y) = trans(p.XL,p.YL)
      p.XL = X
      p.YL = Y
    end
    CSV.write(forinv_data["points"],points)
    rm(forinv_dldirpath,recursive=true)
    return forinv_vars
end
