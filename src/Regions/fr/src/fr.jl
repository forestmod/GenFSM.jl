
module Res_fr
import ..GenFSM.Res

import Downloads 
import ArchGDAL, Rasters, ZipFile, DataStructures # , FTPClient
import Geomorphometry # for slope and aspect 
import Shapefile
import CSV, DataFrames
import Proj # to convert the (X,Y) coordinates of the inventory points
import GeoDataFrames
import Logging

include(joinpath("Res","Utils.jl"))
include(joinpath("Res","Get_data.jl"))
include(joinpath("Res","Init.jl"))

"""
    Res.init!!(::Val{:fr},pixels,settings,overal_region_mask)

Wrapper function to call the initialization function for the French region.

We are in the middle of the following function chain:

Res.init!!() --> init!!(Val(reg_xx)),...) --> Res_xx._init!!()

Note that this function is injected to the `Res` module, but it is defined in the `Res_fr` module (file `fr.jl``)

"""
function Res.init!!(::Val{:fr},pixels,settings,overal_region_mask)
    Res_fr._init!!(pixels,settings,overal_region_mask)
end


end # End module Res_fr


# ------------------------------------------------------------------------------
module Mkt_fr

end

