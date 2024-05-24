module GenFSM

import GenFSM_resources as RES
import GenFSM_simrep as REP

export run

greet() = print("Hello World!")

function run(project="default",scenario="default")
    settings = REP.load_settings(project,scenario)
    years    = settings["years"]
    dc       = settings["dc"]
    ft       = load the data
    region   = settings["simulation_region"]
    init_regions = settings["resource_init_regions"]
    pixels   = RES.init(ft,dc,years,region,init_regions) # return the initialized pixels
end


end # module GenFSM
