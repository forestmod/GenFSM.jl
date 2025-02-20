
export Pixel, make_raster_mask, make_pixels

"""

    $(TYPEDEF)

Pixel used in the Res module



# Fields:
$(TYPEDFIELDS)


"""
mutable struct Pixel
    "A (x,y) tuple. Note that x is longitude and y is latitude (from top left)"
    id::Tuple{Int64,Int64}
    # States
    "Forest volumes, Mm^3, by ft,dc"
    fvol::Array{Float64,2}
    "Forest Area, ha, by ft,dc"
    farea::Array{Float64,2}
    "Availability coefficient, ∈ [0,1] by ft"
    ac::Array{Float64,1} 
    """
    Fixed pixel characteristics

    Encoded characteristics of the pixel (soil, position, geomorphology...).
    These are assumed not to change in the timeframe of the simulation, even if we change FT that may have an effect on soil
    """
    fixed_state::Array{Float64,1}
    """
    Time variable characteristics
    Encoded characteristics of the pixel that by scenario and/or time depends (climate...)
    """
    var_state::Array{Float64,1}
    # Add here some sort of soil degradation that effect depends on the harvesting multiplied by a coefficient of "good vs bad" harvesting. Good harvesting is more expensive but lead to less degradation. 
    """
    Initialized by at least one module
    TODO: maybe unneeded just set initialization_model to 0
    """
    initialized::Bool
    "Initialization model id. It defines the ML model for growth and mortality"
    initialization_model::Int64
    # Forest dynamic parameters
    #mortality_rate::Array{Float64,3} # by ft,dc,t
    #tp::Array{Float64,3} # time of passage to higher DC, by ft,dc,t
end

"""
    make_raster_mask

Create a raster mask of the overall region under study

Note the different layout between the raster geographial representation and the underlying matrix data:

```
Matrix:      Raster:
1 5          1 2 3 4
2 6    ==>   5 6 7 8
3 7
4 8
```

So, when you take raster[2,1] you are taking the raster x=2 (from left), y=1(from top), that is 2 in this example.
"""

function make_raster_mask(region)
    cres_epsg_id = region["cres_epsg_id"]
    x_lb = region["x_lb"]
    x_ub = region["x_ub"]
    y_lb = region["y_lb"]
    y_ub = region["y_ub"]
    xres = region["xres"]
    yres = region["yres"]
    lon, lat = Rasters.X(x_lb:xres:x_ub), Rasters.Y(y_lb:yres:y_ub)
    #crs = convert(Rasters.WellKnownText, Rasters.EPSG(cres_epsg_id)) # dosn't work now, perhaps it needs Proj.jl
    crs = Rasters.EPSG(cres_epsg_id)
    ras = Rasters.Raster(zeros(lon, lat),crs=crs)
    nlon = length(lon)
    nlat = length(lat)
    v = reshape(1:nlon*nlat,nlon,nlat)
    ras = Rasters.Raster(v,dims=(lon,lat),crs=crs) # todo, temporary to check the mask
    ras = Rasters.reverse(ras;dims=Rasters.Y )
    ras[:,:] = v
    return ras
end

"""
   make_pixels(nft,ndc,nx,ny)

Return an array of uninstantiated pixels.

Note that we follow here the geo standard, so pixels[2,1] is pixel x=2 (from left), y=1 (from nord)

"""
function make_pixels(nft,ndc,nx,ny)
    pixels = Array{Pixel,2}(undef,nx,ny) # rows cols

    for ix in 1:nx, iy in 1:ny
        pixels[ix,iy] = Pixel( 
           (ix,iy), # id
           zeros(nft,ndc),      # fvol
           zeros(nft,ndc),      # farea
           zeros(nft), # ac
           zeros(5), # fixed_state
           zeros(5), # var_state
           false,    # initialized
           0,        # initialization_model
        )
    end
    
    return pixels
end

#=
cres_epsg_id: 3035
x_lb: 2.980269469245921e6 # -6 in EPSG:4326 (i.e. WGS 84 ellipsoid unprojected)
x_ub: 4.252324473944207e6 # 9
y_lb: 2.131472196486537e6 # 41
y_ub: 3.210471068338138e6 # 52
xres: 8000 # meters
yres: 8000 # meters
=#

