#"""
#    ResourceInitFrExt
#
#Extension module for initilising the Resource module for France
#"""
#module ResourceInitFrExt 
#
#using GenFSM, GenFSM_resource_init_fr

import GenFSM_resources_init_fr

function resources_init_fr!(pixels,mask,settings)
    # Do the initialization by calling the init function in the package GenFSM_resources_fr
    return GenFSM_resources_init_fr.init!(pixels,mask,settings)
end


#end # module