

GenFSM.runsim
- ScenarioLoader.load_full_settings
- Res.make_raster_mask
- RES.make_pixels
- RES.init!!
  - RES.init!!(Val(fr),...)
    - Res_fr._init!!()
      - Res_fr.get_data!
        - get_dtm
        - get_soil_data
        - get_clc
        - get_climate_data
        - get_inventory_data
      - Res_fr.prepare_data
        - prepare_ign_data
        - train_autoencode_clim_soil
        - train_growth_model
        - define_state
        