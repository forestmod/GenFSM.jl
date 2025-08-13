
module Res_fr
import ..GenFSM
import ..GenFSM: NONE,LOW,STD,HIGH,FULL
import ..GenFSM.Res

using DocStringExtensions, Dates

import Downloads 
import StatsBase
import DataStructures: OrderedDict
import Pipe: @pipe
import HTTP
import ArchGDAL, Rasters, ZipFile, DataStructures # , FTPClient
import Geomorphometry # for slope and aspect 
import Shapefile
import CSV, DataFrames
import Proj # to convert the (X,Y) coordinates of the inventory points
import GeoDataFrames
import Logging
import BetaML
import JLD2

import PythonCall
import CondaPkg
const Chelsa      = PythonCall.pyimport("chelsa_cmip6")
const GetClim     = PythonCall.pyimport("chelsa_cmip6.GetClim")



include(joinpath("Res","Utils.jl"))
include(joinpath("Res","Get_data.jl"))
include(joinpath("Res","Prepare_data.jl"))
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
using DocStringExtensions

end

