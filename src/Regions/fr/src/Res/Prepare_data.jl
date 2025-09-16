# Prepare the data for the French Rres_fr module.

function prepare_data!(settings, mask)
    settings["verbosity"] >= GenFSM.LOW && @info("Preparing res data for the French region.")
    # This function prepares the data for the French region.
    # It is called after the data has been downloaded and layers saved as tiff.
    
    ign_state, ign_growth = prepare_ign_data(settings, mask)
    # the ae_clim and ae_soil models could be used as sort of pretrain and then chained vertically together to form the final growth/mortality model instead of using only the reduced form, but it would takes a lot of computational power in prediction, that is not suitable 
    xclimh_reduced, scaler_clim_m, ae_clim_m  = train_autoencode_clim(settings, mask, ign_growth)
    xclimf_reduced   = GenFSM.Res_fr.predict_autoencoder_clim(settings, mask, scaler_clim_m, ae_clim_m)
    scaler_clim_m, ae_clim_m = nothing, nothing # removing the climatic ae models from memory, as they are not needed anymore
    # Filtering out x,y that are not in the mask (yes, it could have been done earlier..) NO LONGER NEEDEED
    #xclimh_reduced = xclimh_reduced[getindex.(Ref(mask),xclimh_reduced.C,xclimh_reduced.R) .== 1,:]
    #xclimf_reduced = xclimf_reduced[getindex.(Ref(mask),xclimf_reduced.C,xclimf_reduced.R) .== 1,:]
    xfixedpx_reduced       = trainpredict_autoencode_fixedpxdata(settings, mask) # dtm and soil

    sm,mgrs         = GenFSM.Res_fr.train_growth_model(settings, ign_growth, xclimh_reduced, xfixedpx_reduced)
    #res_mortality_m = train_mortality_model(settings, ign_state, ign_growth, xclimh_reduced, xfixedpx_reduced)
    #define_state(settings,mask)
  
    settings["verbosity"] >= STD && @info("DONE preparing data for French region.")
end

"""

    prepare_ign_data(settings)

This function prepares the IGN data for the French region, namely the `dfgrowth` and `dfstate` DataFrames.

**dfgrowth**: 

df: growing_points
filtered by: 
- points with 2 visits
- no interventions/recent cuts
- homogeneous vertical structure

idp: id point
x: x coord
y: y coord
xid: x index in the mask
yid: y index in the mask
y1: year of the first visit
y2: year of the second visit
vHa1: volume per hectare (alive) at the first visit
vHa2: volume per hectare (alive) at the second visit
vHaD: volume of dead trees per hectare in between the two visits
d1: avg diameter at the first visit
d2: avg diameter at the second visit
spg: species group BR,MIX,CON
gov: structure HF,MIX,COP


**dfstate**:

df: starting_points
filtered by: 
- y: latest 3 years

idp: id point
x: x coord
y: y coord
xid: x index in the mask
yid: y index in the mask
year: year
vHa: volume per hectare
d: avg diameter
spg: species group BR,MIX,CON
gov: government HF,MIX,COP

"""
function prepare_ign_data(settings, mask )
 
    verbosity = settings["verbosity"]

    # This function prepares the IGN data for the French region
    verbosity  >= GenFSM.STD && @info("- preparing IGN data for growth/mortality model training and set initial forest resource condition")

    force_other = settings["res"]["fr"]["force_other"]
    force_ml_train = settings["res"]["fr"]["force_ml_train"]
    basefolder  = joinpath(settings["res"]["fr"]["cache_path"],"forinv")
    hyears    = settings["res"]["fr"]["data_sources"]["clim"]["hist_years"]
    ftypes    = settings["res"]["fr"]["ftypes"]

    # if dfgrowth exists and no "ing_points""in force_other, then returning the saved dfgrowth, dfstate
    if (! ("ing_points" in force_other)) && (isfile(joinpath(basefolder,"dfgrowth.csv")) && isfile(joinpath(basefolder,"dfstate.csv"))  )
        dfgrowth = CSV.File(joinpath(basefolder,"dfgrowth.csv")) |> DataFrames.DataFrame
        dfstate = CSV.File(joinpath(basefolder,"dfstate.csv")) |> DataFrames.DataFrame
        return dfstate,dfgrowth 
    end

    tobs  = CSV.File(joinpath(basefolder,"ARBRE.csv")) |> DataFrames.DataFrame
    pobs  = CSV.File(joinpath(basefolder,"PLACETTE.csv")) |> DataFrames.DataFrame
    tpobs = DataFrames.innerjoin(tobs,pobs,on=["IDP","CAMPAGNE"])
    species = CSV.File(joinpath(basefolder,"espar-cdref13.csv")) |> DataFrames.DataFrame
    genders_broadleaves = settings["res"]["fr"]["ign_data"]["genders_broadleaves"]
    genders_coniferous = settings["res"]["fr"]["ign_data"]["genders_coniferous"]
    nsample_v_imputation = settings["res"]["fr"]["ign_data"]["nsample_v_imputation"]
    species.gender = [split(lat_name," ")[1] for lat_name in species.lib_cdref]
    species.group  = [(gender in genders_broadleaves) ? "br" : ((gender in genders_coniferous) ? "con" : missing ) for gender in species.gender] 
    species.spcode = [length(espar) == 1 ? "0$espar" : espar for espar in species[:,"// espar"]] # correction because the use stuff like "2" in the classification, but "02" in the data... arghhh

    tpobs.D = tpobs.C13 ./ pi
    # Form factors
    tpobs.SQFF = tpobs.V ./ tpobs.D .^ 2
    tpobs.CUBFF = tpobs.V ./ tpobs.D .^ 3
    # Vegetation state
    tpobs.VSTATE = [ismissing(t.VEGET5) ? t.VEGET : t.VEGET5  for t in eachrow(tpobs)] 

    # Species group
    map_espar = Dict(species.spcode .=> species.group)
    # Manually add codes that are not in the species list....
    map_espar["332G"] = "br"  # populus
    map_espar["29AF"] = "br"  # other broadleaves
    map_espar["68CE"] = "con" # other coniferous 
    map_espar["25E3"] = "br"  # some salix
    map_espar["25E5"] = "br"  # some salix
    tpobs.SPGR = [ismissing(t.ESPAR) ? missing : map_espar[t.ESPAR] for t in eachrow(tpobs)]

    pids  = unique(pobs.IDP)
    years = unique(pobs.CAMPAGNE)

    map_structure_sfo   = Dict(0 => "", 1 => "HF", 2 => "", 3 => "MIX", 4 => "COP")  # skipping irregular HF (2) and no structure (0)
    map_structure_sver  = Dict("0" => "", "X" => "", "2" => "HF", "3" => "COP", "4" => "", "5" => "MIX", "6" => "HF") # skipping no structure/cutted (0,X) and irregular HF (4)
    map_recent_cuts = Dict(0 => false, 1 => true, 2 => true)

    Xdim = Rasters.dims(mask, Rasters.X)    # a DimensionalData.Sampled{<:Real} or similar
    Ydim = Rasters.dims(mask, Rasters.Y)  
    xcoords = collect(Xdim)
    ycoords = collect(Ydim)
    y_l = length(ycoords)
    y_rev = reverse(ycoords)

    function computeW(c13)
        ismissing(c13) && return missing
        if c13 < 0.705
            return 10000/(6^2*pi)
        elseif c13 < 1.175
            return 10000/(9^2*pi)
        else 
            return 10000/(15^2*pi)
        end
    end    



    # Computing "generic" and sqare and cubic form factors for trees that doesn't have an entry on the first visit, by (1) species, department and diameter class, (2) species and diameter class, (3) species only
    dclasses = [0:10:130;100000]
    tpobs.DC = [ismissing(d) ? missing : (findfirst(x -> x >= d*100, dclasses) -1) for d in tpobs.D]

    ff_sp_dep_dc =DataFrames.combine(DataFrames.groupby(tpobs, [:ESPAR,:SER,:DC,])) do subdf
        # subdf is a DataFrame with the same columns as tpobs
        # we compute the mean of V for each group
        size(subdf,1) == 0 && return (n = 0, sqff = missing, cubff = missing)
        all(ismissing.(subdf.V)) && return (n=0, sqff = missing, cubff = missing)
        (n = size(subdf,1), sqff  = StatsBase.median(skipmissing(subdf.SQFF)),
        cubff = StatsBase.mean(skipmissing(subdf.CUBFF)))
    end
    ff_sp_dep_dc =  GenFSM.to_dict(ff_sp_dep_dc, [:ESPAR,:SER,:DC], [:n,:sqff,:cubff])

    ff_sp_dc =DataFrames.combine(DataFrames.groupby(tpobs, [:ESPAR,:DC,])) do subdf
        # subdf is a DataFrame with the same columns as tpobs
        # we compute the mean of V for each group
        size(subdf,1) == 0 && return (n = 0, sqff = missing, cubff = missing)
        all(ismissing.(subdf.V)) && return (n=0, sqff = missing, cubff = missing)
        (n = size(subdf,1), sqff  = StatsBase.median(skipmissing(subdf.SQFF)),
        cubff = StatsBase.mean(skipmissing(subdf.CUBFF)))
    end
    ff_sp_dc =  GenFSM.to_dict(ff_sp_dc, [:ESPAR,:DC], [:n,:sqff,:cubff])

    ff_sp =DataFrames.combine(DataFrames.groupby(tpobs, [:ESPAR])) do subdf
        # subdf is a DataFrame with the same columns as tpobs
        # we compute the mean of V for each group
        size(subdf,1) == 0 && return (n = 0, sqff = missing, cubff = missing)
        all(ismissing.(subdf.V)) && return (n=0, sqff = missing, cubff = missing)
        (n = size(subdf,1), sqff  = StatsBase.median(skipmissing(subdf.SQFF)),
        cubff = StatsBase.mean(skipmissing(subdf.CUBFF)))
    end
    ff_sp =  GenFSM.to_dict(ff_sp, [:ESPAR], [:n,:sqff,:cubff])


    ff_dep_dc =DataFrames.combine(DataFrames.groupby(tpobs, [:SER,:DC,])) do subdf
        # subdf is a DataFrame with the same columns as tpobs
        # we compute the mean of V for each group
        size(subdf,1) == 0 && return (n = 0, sqff = missing, cubff = missing)
        all(ismissing.(subdf.V)) && return (n=0, sqff = missing, cubff = missing)
        (n = size(subdf,1), sqff  = StatsBase.mean(skipmissing(subdf.SQFF)),
        cubff = StatsBase.mean(skipmissing(subdf.CUBFF)))
    end
    ff_dep_dc =  GenFSM.to_dict(ff_dep_dc, [:SER,:DC], [:n,:sqff,:cubff])

    missingcheck(res) = ismissing(res) || ismissing(res[2])

    function get_generic_ff(sp,dep,dc)
        # returns the generic form factor for a given species, department and diameter class
        # if not found, returns missing
        (n,sqff,cubff) = (0,0.0,0.0)
        if ismissing(sp) 
            (n,sqff,cubff) = (ismissing(get(ff_dep_dc,(dep,dc),missing) )) ? (missing,missing,missing) : get(ff_dep_dc,(dep,dc),missing)
        elseif ! missingcheck(get(ff_sp_dep_dc,(sp, dep,dc),missing) )
            (n,sqff,cubff) =  get(ff_sp_dep_dc,(sp, dep,dc),missing)
        elseif ! missingcheck(get(ff_sp_dc,(sp,dc),missing) )
            (n,sqff,cubff) =  get(ff_sp_dc,(sp,dc),missing)
        elseif ! missingcheck(get(ff_sp,(sp),missing) )
            (n,sqff,cubff) =  get(ff_sp,(sp),missing)
        else
            (n,sqff,cubff) =  missing,missing, missing
        end
        return (!ismissing(sqff) && sqff < 20) ? (n,sqff,cubff) : (missing, missing, missing)
    end


    # Computing, after imputation, the vHa for each tree
    #tpobs.vHa  = vHaContribution.(tpobs.V,tpobs.C13)

    # Selecting only points with 2 visits
    cmap = StatsBase.countmap(pobs.IDP)
    pids2 = filter(x -> cmap[x] == 2, pids)
    #pobs2  = pobs[in.(pobs.IDP, Ref(pids2)), :]

    # Step 1: computing dfgrowth data
    @info("Starting building dfgrowth DataFrame...")
    skipped = Dict("harvesting_detected" => 0,
                "out_of_mask" => 0,
                "out_of_time_range" => 0,
                "unsimulated_ftype" => 0,
                "no_vertical_structure" => 0,
                "missing_data" => 0,
                "unreliable_data1" => 0,
                "unreliable_data2" => 0,
                "unreliable_data3" => 0,
                "empty_plots" => 0,
    )

    dfgrowth = DataFrames.DataFrame(idp = Int64[], x = Float64[], y = Float64[], xid = Int64[], yid = Int64[], y1 = Int64[], y2 = Int64[], vHa1 = Float64[], vHa2 = Float64[], vHaD = Float64[], d1 = Float64[], d2 = Float64[], spg = String[], gov = String[])
    for (i,p) in enumerate(pids2)
        #i == 10000 && break
        #i = 500098
        p = pids2[i]
        #println("p: $p")
        #p = 774806 # 507076 #500222
        #i= 0
        i%1000 == 0 && println("p: $p - Processed points: $i/$(length(pids2)) - Added points $(size(dfgrowth,1)) - skipped: $skipped")
        #p = pids2[end]
        #p = 146283 
        visit1   = pobs[pobs.IDP .== p .&& pobs.VISITE .== 1, :][1,:]
        visit2   = pobs[pobs.IDP .== p .&& pobs.VISITE .== 2, :][1,:]

        # Filters 1: early filters 
        xid = searchsortedlast.(Ref(xcoords), visit1.XL) 
        yid = y_l .- searchsortedlast.(Ref(y_rev), visit1.YL)
        getindex(mask,xid,yid) == 1 || (skipped["out_of_mask"] += 1; continue ) # skip points outside the study area
        # Vertical structure
        vs= ""
        if !ismissing(visit1.SFO)
            vs = map_structure_sfo[visit1.SFO]
        elseif !ismissing(visit1.SVER) 
            vs = map_structure_sver[visit1.SVER]
        end
        vs !== "" || (skipped["no_vertical_structure"] +=1; continue) # skip points with no vertical structure
        (in(visit1.CAMPAGNE,hyears) && in(visit2.CAMPAGNE,hyears)) || (skipped["out_of_time_range"] += 1; continue ) # skip points outside the time range

        treefields = ["VISITE", "A", "VSTATE", "ESPAR", "SPGR", "SER", "C13", "D", "DC", "SQFF", "CUBFF", "V", "W"]
        mytrees  = tpobs[tpobs.IDP .== p, treefields]
        mytrees1 = @view mytrees[mytrees.VISITE .==1, :]
        mytrees2 = @view mytrees[mytrees.VISITE .==2, :]
        
        # Filters 2: mid-time filters
        any(in.(mytrees2.VSTATE, Ref(["6","7"]))) && (skipped["harvesting_detected"] += 1; continue ) # skip points with harvesting detected
        size(mytrees2,1) < 2  && (skipped["empty_plots"] += 1; continue )   # skip points with < 2 (dead or alive) trees on the second visit
        (length(collect(skipmissing(mytrees.SQFF))) > 0 &&  maximum(collect(skipmissing(mytrees.SQFF))) < 20) || (skipped["unreliable_data1"] += 1; continue )   # skip points with unreliable ratio V/D
        (length(collect(skipmissing(mytrees.SQFF))) > 0 && minimum(collect(skipmissing(mytrees.SQFF))) > 0.1) || (skipped["unreliable_data2"] += 1; continue )  # skip points with unreliable ratio V/D

        mytrees21 = DataFrames.leftjoin(mytrees2, mytrees1, on = :A, makeunique=true)

        mytrees21_alive  = @view mytrees21[mytrees21.VSTATE .== "0", :]
        mytrees21_alive1 = @view mytrees21[.!ismissing.(mytrees21.VSTATE_1) .&& mytrees21.VSTATE_1 .== "0", :]
        mytrees21_dead   = @view mytrees21[( (.! in.(mytrees21.VSTATE, Ref(["0","6","7"]))) # new dead only (not dead in the first visit)
                                        .&& 
                                    (ismissing.(mytrees21.VISITE_1) 
                                        .|| 
                                        (.! ismissing.(mytrees21.VSTATE_1)
                                    .&& mytrees21.VSTATE_1 .== "0"))),:]
        # Assigning to dead trees the diameter of first visit if missed
        [t.D .== t.D_1 for t in eachrow(mytrees21_dead) if ismissing(t.D) && !ismissing(t.VISITE_1) && !ismissing(t.D_1)]


        # Filter 3: final filter
        # All alive2 trees and new-entry dead trees must have diameter
        (sum(ismissing.(mytrees21_alive.D)) > 0 || sum(ismissing.(mytrees21_dead.D)) > 0 ) && (skipped["missing_data"] += 1; continue; )
        sum(ismissing.(mytrees21_alive1.D_1)) > 0  && (skipped["missing_data"] += 1; continue; )
        all(ismissing.(mytrees21_alive.VISITE_1) .|| ( mytrees21_alive.D .>= mytrees21_alive.D_1 .&&  (mytrees21_alive.D .<= 0.15 .|| mytrees21_alive.D .<= 1.5 * mytrees21_alive.D_1)          ) ) || (skipped["unreliable_data3"] += 1; continue  ) 
        all(ismissing.(mytrees21_dead.VISITE_1) .|| ( mytrees21_dead.D .>= mytrees21_dead.D_1 .&&  (mytrees21_dead.D .<= 0.15 .|| mytrees21_dead.D .<= 1.5 * mytrees21_dead.D_1)          ) ) || (skipped["unreliable_data3"] += 1; continue)   
        sum(mytrees21_alive1.D_1) > 0.05 || (skipped["empty_plots"] += 1; continue)  # otherwise difficult to make exponential growth

        # Assigning W factor
        [t.W = t.W_1 for t in eachrow(mytrees21) if ismissing(t.W) && !ismissing(t.VISITE_1) && !ismissing(t.W_1)]
        [t.W = computeW(t.C13) for t in eachrow(mytrees21) if ismissing(t.W) && (ismissing(t.VISITE_1) || ismissing(t.W_1))]

        # Assigning species if missing but species not missing in first visit
        [t.ESPAR = t.ESPAR_1 for t in eachrow(mytrees21) if ismissing(t.ESPAR) && !ismissing(t.VISITE_1) && !ismissing(t.ESPAR_1)]

        # computing species group
        [t.SPGR = ismissing(t.ESPAR) ? missing : map_espar[t.ESPAR] for t in eachrow(mytrees21)]


        # Assigning volumes, vHA and species groups
        [t.V = 0.8 * t.D^3* t.CUBFF_1 + 0.2 * t.D^2* t.SQFF_1 for t in eachrow(mytrees21) if ismissing(t.V) && !ismissing(t.VISITE_1) && !ismissing(t.D_1) && !ismissing(t.V_1)]
        [(ff = get_generic_ff(t.ESPAR,t.SER,t.DC); t.V = 0.8 * t.D^3* ff[3] + 0.2 * t.D^2* ff[2])  for t in eachrow(mytrees21) if ismissing(t.V) && (ismissing(t.VISITE_1) || ismissing(t.D_1) || ismissing(t.V_1))]

        mytrees21.vHa = mytrees21.V .* mytrees21.W
        mytrees21.vHa_1 = mytrees21.V_1 .* mytrees21.W_1

        # Defining the species group
    
        vHaBr  = sum( mytrees21_alive[mytrees21_alive.SPGR .== "br", "vHa"])
        vHaCon = sum( mytrees21_alive[mytrees21_alive.SPGR .== "con", "vHa"])
        vHaTot = vHaBr .+ vHaCon
        vHaBrRatio = vHaBr ./ vHaTot
        if vHaBrRatio < 0.1 
            spgr = "con"
        elseif vHaBrRatio < 0.9
            spgr = "mixsp"
        else
            spgr = "br"
        end
        vHa1 = sum(mytrees21_alive1.vHa_1)
        vHa2 = sum(mytrees21_alive.vHa)
        vHaD = sum(mytrees21_dead.vHa)
        d1 = (size(mytrees21_alive1,1) == 0) ? 0.0 : StatsBase.mean(mytrees21_alive1.C13_1)  # some plots have no trees on the first visit, so no diameter. This is ok when the trees are too small to be measured
        d2 = (size(mytrees21_alive,1) == 0) ? d1 : StatsBase.mean(mytrees21_alive.C13)
        # Checking for no missing values...
        ismissing(visit1.XL) && error("XL missing for point $p")
        ismissing(visit1.YL) && error("YL missing for point $p")
        ismissing(visit1.CAMPAGNE) && error("CAMPAGNE1 missing for point $p")
        ismissing(visit2.CAMPAGNE) && error("CAMPAGNE2 missing for point $p")
        ismissing(vHa1) && error("vHa1 missing for point $p")
        ismissing(vHa2) && error("vHa2 missing for point $p")
        ismissing(vHaD) && error("vHaD missing for point $p")
        ismissing(d1) && error("d1 missing for point $p")
        ismissing(d2) && error("d2 missing for point $p")
        ismissing(spgr) && error("spgr missing for point $p")
        ismissing(vs) && error("vs missing for point $p")
        (vHa2+vHaD) <  (vHa1 - 0.00000001) && error("vHa2 + vHaD < vHa1 for point $p")
        # Now I have everything I need to create the row for the growing_points DataFrame

        push!(dfgrowth, [p, visit1.XL, visit1.YL, xid, yid, visit1.CAMPAGNE, visit2.CAMPAGNE, vHa1, vHa2, vHaD, d1, d2, spgr, vs])
    end


    unmatch = size(pids2,1)- sum(skipped[k] for k in keys(skipped)) - size(dfgrowth,1) 
    verbosity >= GenFSM.HIGH && println("Unmatch in dfgrowth df creation (should be zero): $(unmatch)")
    verbosity >= GenFSM.HIGH && println(skipped)

    dfgrowth.dV     = (dfgrowth.vHa2+dfgrowth.vHaD) .- dfgrowth.vHa1
    dfgrowth.dV2    = dfgrowth.vHa2 .- dfgrowth.vHa1
    dfgrowth.mShare = dfgrowth.vHaD ./ (dfgrowth.vHa2 .+ dfgrowth.vHaD) # share of dead volume in the total volume change
    dfgrowth.ftype = ["$(r.spg)_$(r.gov)" for r in eachrow(dfgrowth)]

    # Filter the volume growth df by observed historical weather, forest type simulmated and spatially within the mask grid
    dfgrowth = dfgrowth[in.(dfgrowth.y1,Ref(hyears)) .&& in.(dfgrowth.y2,Ref(hyears)) .&& in.(dfgrowth.ftype,Ref(ftypes)),:]
    #dfgrowth.xid = searchsortedlast.(Ref(xcoords), dfgrowth.x) 
    #dfgrowth.yid = y_l .- searchsortedlast.(Ref(y_rev), dfgrowth.y)
    #dfgrowth = dfgrowth[getindex.(Ref(mask),dfgrowth.xid,dfgrowth.yid) .== 1,:]


    CSV.write(joinpath(basefolder,"dfgrowth.csv"),dfgrowth)


    @info("Done building dfgrowth DataFrame.")

    # Step 2: building dfstate DataFrame
    @info("Starting building dfstate DataFrame...")
    pid_latest = pobs[ in.(pobs.CAMPAGNE, Ref(maximum(pobs.CAMPAGNE)-2 : maximum(pobs.CAMPAGNE))), "IDP"]
    #443593 in first visit and 218916 in second visit

    counts = fill(0,4)
    dfstate = DataFrames.DataFrame(idp = Int64[], x = Float64[], y = Float64[], xid=Float64[], yid=Float64[],year=Int64[], vHa = Float64[], d = Float64[], spg = String[], gov = String[])
    for (i,p) in enumerate(pid_latest)
        counts[1] += 1
        i%1000 == 0 && println("p: $p np: $(length(pid_latest)) counts: $counts")
        #p = pid_latest[end]
        #p = 146283 
        visit = pobs[pobs.IDP .== p, :][1,:]

        xid = searchsortedlast.(Ref(xcoords), visit.XL) 
        yid = y_l .- searchsortedlast.(Ref(y_rev), visit.YL)
        getindex(mask,xid,yid) == 1 || continue # skip points outside the study area

        in.(visit.CSA, Ref(["1","3","5","2"])) || continue # skip points with no forest structure (open, closed forest or popplars)
        counts[2] += 1

        mytrees = tpobs[tpobs.IDP .== p, :]
        mytrees.vegstate = [ismissing(t.VEGET5) ? t.VEGET : t.VEGET5  for t in eachrow(mytrees)] # some second visit trees have veg state in VEGET instead of VEGET5 :-()
        mytrees_alive  = mytrees[mytrees.vegstate .== "0", :]

        # estimating V, but only if the FF of the tree has already been computed on a previous occasion
        for (j,t) in enumerate(eachrow(mytrees_alive))
            if ismissing(t.V) && t.VISITE == 2
                t1s = tpobs[tpobs.IDP .== t.IDP .&& tpobs.A .== t.A .&& tpobs.VISITE .== 1, :]
                if size(t1s,1) == 1 && !ismissing(t1s[1,"V"]) # there is a match on first visit and the match has a volume
                    t1 = t1s[1,:]
                    sqFF, cubFF = (t1.SQFF, t1.CUBFF)
                    # if also the diameter is missing, we give it the diameter of the previous 5 years
                    if ismissing(t.D)
                        t.D = t1.D
                        t.C13 = t1.C13
                    end
                    t.V = 0.8 * t.D^3* cubFF + 0.2 * t.D^2* sqFF
                end
            end
            if ismissing(t.W) && t.VISITE == 2
                t1s = tpobs[tpobs.IDP .== t.IDP .&& tpobs.A .== t.A .&& tpobs.VISITE .== 1, :]
                if size(t1s,1) == 1 && !ismissing(t1s[1,"W"]) # there is a match on first visit and the match has a volume
                    t1 = t1s[1,:]
                    # if also the diameter is missing, we give it the diameter of the previous 5 years
                    if ismissing(t.D)
                        t.D = t1.D
                        t.C13 = t1.C13
                    end
                    t.W = t1.W
                end
            elseif ismissing(t.W) && t.VISITE == 1
                t.W = computeW(t.C13) # computing W parameter
            end

            if t.VISITE == 2
                t1s = tpobs[tpobs.IDP .== t.IDP .&& tpobs.A .== t.A .&& tpobs.VISITE .== 1, :]
                # setting C13, D and V as minima as the value of first visit to reduce errors (e.g. point 507076, tree 11)
                if size(t1s,1) > 0
                    t1   = t1s[1,:]
                    t.C13 = max(t.C13, t1.C13)
                    t.D   = max(t.D, t1.D)
                    t.V   = max(t.V, t1.V)
                end
            end
        end
        #mytrees_alive.vHa = vHaContribution.(mytrees_alive.V, mytrees_alive.C13) # computing vHa for alive trees on second visit
        mytrees_alive.vHa = mytrees_alive.V .* mytrees_alive.W # computing vHa for alive trees on second visit 

        vHa = sum(mytrees_alive.vHa)

        ismissing(vHa) && continue # skip points with no volume (or diameterinfo on alive trees
        counts[3] += 1


        # Defining the species group
        mytrees_alive.spgr = [ismissing(t.ESPAR) ? "" : map_espar[t.ESPAR] for t in eachrow(mytrees_alive)]
        vHaBr  = sum( mytrees_alive[mytrees_alive.spgr .== "br", "vHa"])
        vHaCon = sum( mytrees_alive[mytrees_alive.spgr .== "con", "vHa"])
        vHaTot = vHaBr .+ vHaCon
        vHaBrRatio = vHaBr ./ vHaTot
        if vHaBrRatio < 0.1 
            spgr = "con"
        elseif vHaBrRatio < 0.9
            spgr = "mixsp"
        else
            spgr = "br"
        end

        # Vertical structure
        vs= ""
        if !ismissing(visit.SFO)
            vs = map_structure_sfo[visit.SFO]
        elseif !ismissing(visit.SVER) 
            vs = map_structure_sver[visit.SVER]
        end
        vs !== "" || continue # skip points with no vertical structure
        counts[4] += 1

        d = (size(mytrees_alive,1) == 0) ? 0.0 : StatsBase.mean(mytrees_alive.C13)  # some plots have no trees so no diameter. This is ok when the trees are too small to be measured or has just been cutted

        # Checking for no missing values...
        ismissing(visit.XL) && error("XL missing for point $p")
        ismissing(visit.YL) && error("YL missing for point $p")
        ismissing(visit.CAMPAGNE) && error("CAMPAGNE missing for point $p")
        ismissing(vHa) && error("vHa missing for point $p")
        ismissing(d) && error("d missing for point $p")
        ismissing(spgr) && error("spgr missing for point $p")
        ismissing(vs) && error("vs missing for point $p")
        # Now I have everything I need to create the row for the growing_points DataFrame
        push!(dfstate, [p, visit.XL, visit.YL, xid, yid, visit.CAMPAGNE, vHa, d, spgr, vs])
    end
    @info("Done building dfgrowth DataFrame.")

    #dfstate.xid = searchsortedlast.(Ref(xcoords), dfstate.x) 
    #dfstate.yid = y_l .- searchsortedlast.(Ref(y_rev), dfstate.y)
    #dfstate = dfstate[getindex.(Ref(mask),dfstate.xid,dfstate.yid) .== 1,:]

    # Saving the growing points data
    CSV.write(joinpath(basefolder,"dfstate.csv"),dfstate)

    return dfstate, dfgrowth 


end


function train_autoencode_clim(settings, mask, ign_growth)
    # This function trains the autoencoder for climatic data
    # It returns the autoencoded climatic data.

    settings["verbosity"] >= STD && @info("- autoencoding climatic data")

    force_other    = settings["res"]["fr"]["force_other"]
    force_ml_train = settings["res"]["fr"]["force_ml_train"]
    basefolder     = joinpath(settings["res"]["fr"]["cache_path"],"clim_historical")

    # if all needed files exists and no "xclim" in force_other, then returning the saved data
    if (! ("ae_climh" in force_ml_train)) && (isfile(joinpath(basefolder,"xclimh_reduced.csv.gz")) && isfile(joinpath(basefolder,"scaler_clim_m.jld2"))  && isfile(joinpath(basefolder,"ae_clim_m.jld2")) )
        xclimh_reduced = CSV.File(joinpath(basefolder,"xclimh_reduced.csv.gz")) |> DataFrames.DataFrame
        scaler_clim_m  = BetaML.model_load(joinpath(basefolder,"scaler_clim_m.jld2"),"ms")
        ae_clim_m      = BetaML.model_load(joinpath(basefolder,"ae_clim_m.jld2"),"ma")
        return xclimh_reduced, scaler_clim_m, ae_clim_m 
    end

    datafiles = settings["res"]["fr"]["input_rasters"]["clim"]["historical"]
    vars      = settings["res"]["fr"]["data_sources"]["clim"]["vars"]
    hyears    = settings["res"]["fr"]["data_sources"]["clim"]["hist_years"]
    ign_minyear = minimum(ign_growth.y1)
    # we select only the years that are greater or equal to the minimum year of the ign data
    hyears    = hyears[hyears .>= ign_minyear]
    ae_nyears = settings["res"]["fr"]["data_sources"]["clim"]["ae_nyears"]
    ae_nsample = settings["res"]["fr"]["data_sources"]["clim"]["ae_nsample"]
    ae_base_nepochs = settings["res"]["fr"]["data_sources"]["clim"]["ae_base_nepochs"]
    ae_max_ntrains  = settings["res"]["fr"]["data_sources"]["clim"]["ae_max_ntrains"]
    ae_hidden_layer_size = settings["res"]["fr"]["data_sources"]["clim"]["ae_hidden_layer_size"]
    ae_encoded_size = settings["res"]["fr"]["data_sources"]["clim"]["ae_encoded_size"]
    verbosity = settings["verbosity"]
    nC,nR     = size(mask)
    #nxclimh   = nC*nR*(length(hyears)-ae_nyears+1)
    nxclimh   = sum(mask)*(length(hyears)-ae_nyears+1)
  
    if (! ("xclimh" in force_other)) && (isfile(joinpath(basefolder,"xclimh.csv.gz")))
        verbosity >= STD && @info(" -- reading xclimh df from saved CSV file")
        xclimh = CSV.File(joinpath(basefolder,"xclimh.csv.gz")) |> DataFrames.DataFrame
    else
        verbosity >= STD && @info(" -- creating xclimh df from raster files")
        xnames   = collect(keys(datafiles))
        n_xnames  = length(xnames)
        xrasters = OrderedDict{Tuple{String, Int64, Int64},Rasters.Raster}([i => Rasters.Raster(datafiles[i]) |> replace_missing for i in xnames])
        # by time (month, year) var name
        xclimh   = DataFrames.DataFrame("C"=>Array{Int64}(undef,nxclimh), "R"=>Array{Int64}(undef,nxclimh), "Y"=>Array{Int64}(undef,nxclimh), vec(["$(v)_m$(m)_yl$(yl)" => Array{Float64}(undef,nxclimh) for m in 1:12, yl in (ae_nyears-1):-1:0, v in vars])...)

        # Reformatting the data in a large matrix (records [pixels] x variables)
        ridx = 1
        for y in hyears[ae_nyears:end]
          verbosity > HIGH && @info("    -- creating xclimh data for year: $y")
          for c in 1:nC
              for r in 1:nR
                    mask[c,r] == 1 || continue # skip pixels outside the mask
                    xclimh[ridx,[1,2,3]] .= [c,r,y]
                    cidx = 4
                    for v in vars
                        for l in (ae_nyears-1):-1:0
                            for m in 1:12
                                 xclimh[ridx,cidx] = xrasters[(v,m,y-l)][c,r]
                                cidx +=1
                            end
                        end
                    end
                    ridx +=1
              end
            end
        end
        CSV.write(joinpath(basefolder,"xclimh.csv.gz"),xclimh;compress=true)
    end


    srows = StatsBase.sample(axes(xclimh, 1), ae_nsample,replace = false)
    sx    = view(xclimh, srows, :)
    datax = Matrix(sx[:,4:end])
    ((xtrain,xval),) = BetaML.partition([datax],[0.8,0.2])

    ms    = BetaML.Scaler(cache=false)

    BetaML.fit!(ms,xtrain)
    xtrains  = BetaML.predict(ms,xtrain)
    BetaML.model_save(joinpath(basefolder,"scaler_clim_m.jld2");ms)


    ma    = BetaML.AutoEncoder(encoded_size=ae_encoded_size,layers_size=ae_hidden_layer_size,epochs=ae_base_nepochs,verbosity=BetaML.LOW, cache=false)

    last_mse_val = Inf
    xclimh_reduced = nothing

    for n in 1:ae_max_ntrains
        verbosity >= STD && @info("  -- training passage: $n")
        previous_model = deepcopy(ma)
        BetaML.fit!(ma,xtrains)
        xtrainr  = BetaML.predict(ma,xtrains)
        x̂trains  = BetaML.inverse_predict(ma,xtrainr)
        x̂train   = BetaML.inverse_predict(ms,x̂trains)
        mses_train  = [BetaML.mse(xtrain[:,i],x̂train[:,i]) for i in axes(xtrain,2)]
        #means_train = mean(xtrain,dims=1)'
        xvals    = BetaML.predict(ms,xval)
        xvalsr   = BetaML.predict(ma,xvals)
        x̂vals    = BetaML.inverse_predict(ma,xvalsr)
        x̂val     = BetaML.inverse_predict(ms,x̂vals)
        mses_val = [BetaML.mse(xval[:,i],x̂val[:,i]) for i in axes(xval,2)]
        #means_val = mean(xval,dims=1)'
        mse_train_overall = StatsBase.mean(mses_train)
        mse_val_overall   = StatsBase.mean(mses_val)
        verbosity >= STD && @info "  -- $(n): mse_train: $(mse_train_overall)"
        verbosity >= STD && @info "  -- $(n): mse_val: $(mse_val_overall)"
        if mse_val_overall > last_mse_val
            verbosity >= STD && @info "- DONE ae climatic model training."
            verbosity > STD && @info  "  -- detailed mse_train: $(mses_train')"
            verbosity > STD && @info  "  -- detailed mse_val: $(mses_val')"
            ma = previous_model
            BetaML.model_save(joinpath(basefolder,"ae_clim_m.jld2");ma)
            break
        else
            verbosity >= STD && @info "Further training needed..."
            last_mse_val = mse_val_overall
            if n == ae_max_ntrains
                verbosity >= LOW && @warn "Maximum number of training iterations reached, validation still descending. Saving the AE model, but consider further training."
                BetaML.model_save(joinpath(basefolder,"ae_clim_m.jld2");ma)
            end
        end
    end

    xtot = BetaML.predict(ms,Matrix(xclimh[:,4:end]))
    xtot_reduced = BetaML.predict(ma,xtot)
    xclimh_reduced= hcat(xclimh[:,1:3],DataFrames.DataFrame(xtot_reduced,:auto))
    CSV.write(joinpath(basefolder,"xclimh_reduced.csv.gz"),xclimh_reduced;compress=true)

    return (xclimh_reduced,ms,ma)
end


function predict_autoencoder_clim(settings, mask, scaler_clim_m, ae_clim_m)

    settings["verbosity"] >= STD && @info("- predicting autoencoded future climatic data for scenario $(settings["scenario"])")

    force_other    = settings["res"]["fr"]["force_other"]
    basefolder_h     = joinpath(settings["res"]["fr"]["cache_path"],"clim_historical")
    basefolder_f     = joinpath(settings["res"]["fr"]["cache_path"],"clim_future",settings["scenario"])

    # if all needed files exists and no "xclim" in force_other, then returning the saved data
    if (! ("xclimf_reduced" in force_other)) && (isfile(joinpath(basefolder_f,"xclimf_reduced.csv.gz")) )
        xclimf_reduced = CSV.File(joinpath(basefolder_f,"xclimf_reduced.csv.gz")) |> DataFrames.DataFrame
        return xclimf_reduced
    end

    fyears    = settings["res"]["fr"]["data_sources"]["clim"]["fut_years"]
    ae_nyears = settings["res"]["fr"]["data_sources"]["clim"]["ae_nyears"]

    datafiles_h = settings["res"]["fr"]["input_rasters"]["clim"]["historical"]
    datafiles_f = settings["res"]["fr"]["input_rasters"]["clim"]["future"]
    vars      = settings["res"]["fr"]["data_sources"]["clim"]["vars"]
    verbosity = settings["verbosity"]
    nC,nR     = size(mask)
    #nxclimf  = nC*nR*(length(fyears))
    nxclimf   = sum(mask)*(length(fyears))
    ae_encoded_size = settings["res"]["fr"]["data_sources"]["clim"]["ae_encoded_size"]

    xnames_h   = collect(keys(datafiles_h))
    xnames_f   = collect(keys(datafiles_f))
    xrasters_h = OrderedDict{Tuple{String, Int64, Int64},Rasters.Raster}([i => Rasters.Raster(datafiles_h[i]) |> Rasters.replace_missing for i in xnames_h])
    xrasters_f = OrderedDict{Tuple{String, Int64, Int64},Rasters.Raster}([i => Rasters.Raster(datafiles_f[i]) |> Rasters.replace_missing for i in xnames_f])

    xclimf_reduced   = DataFrames.DataFrame("C"=>Array{Int64}(undef,nxclimf), "R"=>Array{Int64}(undef,nxclimf), "Y"=>Array{Int64}(undef,nxclimf), vec(["x$(i)" => Array{Float64}(undef,nxclimf) for i in 1:ae_encoded_size])...)

    ridx = 1
    for y in fyears
    verbosity >= HIGH && @info("    -- creating xclimf_reduced data for year: $y")
    for c in 1:nC
        for r in 1:nR
                mask[c,r] == 1 || continue # skip pixels outside the mask
                cidx = 1
                clim_data_y_px = Array{Float64,1}(undef,ae_nyears*12*length(vars))
                for v in vars
                    for l in (ae_nyears-1):-1:0
                        for m in 1:12
                            if (y-l) < fyears[1]
                            clim_data_y_px[cidx] = xrasters_h[(v,m,y-l)][c,r]
                            else
                            clim_data_y_px[cidx] = xrasters_f[(v,m,y-l)][c,r]
                            end
                            cidx +=1
                        end
                    end
                end
                # ae predict 240 -> 6
                clim_data_red   =  @pipe clim_data_y_px|> transpose |>  BetaML.predict(scaler_clim_m, _) |> BetaML.predict(ae_clim_m, _) |> vec
                xclimf_reduced[ridx,:] .= vcat([c,r,y],clim_data_red)
                ridx +=1
        end
        end
    end
    CSV.write(joinpath(basefolder_f,"xclimf_reduced.csv.gz"),xclimf_reduced;compress=true)
    return xclimf_reduced
end


function trainpredict_autoencode_fixedpxdata(settings, mask)
    # function trainpredict_autoencode_fixedpxdata(settings, mask)
    # Note: we put together here both soil and elevation data.. perhaps it is better to separate them, as dtm is not really a soil variable, but rather a topographic variable, and the ae doesn't work supergood with it.
    verbosity = settings["verbosity"]
    verbosity >= GenFSM.STD && @info("- autoencoding px fixed data")

    force_other    = settings["res"]["fr"]["force_other"]
    basefolder     = joinpath(settings["res"]["fr"]["cache_path"],"pxfixed")

    # if all needed files exists and no "xclim" in force_other, then returning the saved data
    if (! ("xfixedpx_reduced" in force_other)) && isfile(joinpath(basefolder,"xfixedpx_reduced.csv.gz"))
        xfixedpx_reduced = CSV.File(joinpath(basefolder,"xfixedpx_reduced.csv.gz")) |> DataFrames.DataFrame
        return xfixedpx_reduced
    end
    isdir(basefolder) || mkpath(basefolder)

    datafiles_dtm = settings["res"]["fr"]["input_rasters"]["dtm"]
    datafiles_soil = settings["res"]["fr"]["input_rasters"]["soil"]
    datafiles = merge(datafiles_soil,datafiles_dtm) # first soil, so I have the classes not to scalar at the beginning
    ae_maxntrain = settings["res"]["fr"]["px_fixed_data"]["ae_maxntrain"]
    ae_base_nepochs = settings["res"]["fr"]["px_fixed_data"]["ae_base_nepochs"]
    ae_max_ntrains  = settings["res"]["fr"]["px_fixed_data"]["ae_max_ntrains"]
    ae_hidden_layer_size = settings["res"]["fr"]["px_fixed_data"]["ae_hidden_layer_size"]
    ae_encoded_size = settings["res"]["fr"]["px_fixed_data"]["ae_encoded_size"]
    xnames    = collect(keys(datafiles))
    n_xnames  = length(xnames)
    nC,nR = size(mask)
    nx    = sum(mask)
    #nd_dtm  = length(keys(datafiles_dtm))
    #nd_soil = length(keys(datafiles_soil))
    soil_texture_n_classes = settings["res"]["fr"]["data_sources"]["soil_texture_n_classes"]

    verbosity >= GenFSM.STD && @info(" -- creating xfixedpx df from raster files")

    xrasters = OrderedDict{String,Rasters.Raster}([i => Rasters.Raster(datafiles[i]) |> Rasters.replace_missing for i in xnames])
    # by time (month, year) var name
    xfixedpx   = DataFrames.DataFrame("C"=>Array{Int64}(undef,nx), "R"=>Array{Int64}(undef,nx), vec(["$(var)" => Array{Float64}(undef,nx) for var in xnames])...)

    # Reformatting the data in a large matrix (records [pixels] x variables)
    ridx = 1
    (n_problematic_px, n_very_problematic_px) = (0,0)
    for c in 1:nC
        verbosity > GenFSM.HIGH && @info("    -- creating xfixedpx  data for column: $c")
        for r in 1:nR
            mask[c,r] == 1 || continue # skip pixels outside the mask
            xfixedpx[ridx,[1,2]] .= [c,r]
            cidx = 3
            for v in xnames
                if(ismissing(xrasters[v][c,r]))
                    # This happens in a few pixels at the border of the mask. Assigning the average of the variable in the neighmour pixels
                    vals = Float64[]
                    for c2 in max(1,c-2):min(nC,c+2)
                        for r2 in max(1,r-2):min(nR,r+2)
                            if(mask[c2,r2] == 1 && !ismissing(xrasters[v][c2,r2]))
                                push!(vals,xrasters[v][c2,r2])
                            end
                        end
                    end
                    if(length(vals) > 0)
                        n_problematic_px += 1
                        xfixedpx[ridx,cidx] = StatsBase.mean(vals)
                    else
                        n_very_problematic_px += 1
                        xfixedpx[ridx,cidx] = StatsBase.mean(skipmissing(xrasters[v]))
                        verbosity >= GenFSM.HIGH && @warn("Missing value for variable `$(v)` at pixel $(c),$(r). No neighbour pixels found to calculate the average. using the mean of all pixels: $(xfixedpx[ridx,cidx])")
                    end
                else
                    #println("(v,c,r): $(v), $c, $r  - $(xrasters[v][c,r])")
                    xfixedpx[ridx,cidx] = xrasters[v][c,r]
                end
                cidx +=1
            end
            ridx +=1
        end
    end

    npx         = size(xfixedpx,1)
    train_ratio = min(ae_maxntrain/npx, 0.8)
    verbosity >= GenFSM.HIGH && @info("    -- finished creating xfixedpx data. Number of pixels with missing values: $(n_problematic_px), number of pixels with no neighbour pixels: $(n_very_problematic_px), number of total pixels: $(ridx)")

    #=
    srows = StatsBase.sample(axes(xfixedpx, 1), ae_nsample,replace = false)
    sx    = view(xfixedpx, srows, :)
    datax = Matrix(sx[:,3:end]) # NOTE: I skip the first two columns (C,R) !!
    nd = size(datax,2)
    =#


    # we don't use C/R in fixedpx
    ((xtrain,xval),(trainids,valids)) = BetaML.partition([Matrix(xfixedpx[:,3:end]),hcat(1:npx)],[train_ratio,1-train_ratio])

    CSV.write(joinpath(basefolder,"xfixedpx.csv.gz"),xfixedpx;compress=true)
    CSV.write(joinpath(basefolder,"xfixedpx_trainids.csv"),Tables.table(Int.(trainids)))
    CSV.write(joinpath(basefolder,"xfixedpx_valids.csv"),Tables.table(Int.(valids)))

    ms    = BetaML.Scaler(cache=false,skip=1:soil_texture_n_classes) # skip to scale the categorical TextureUSDA category cols

    BetaML.fit!(ms,Matrix(xfixedpx[:,3:end])) # old: BetaML.fit!(ms,xtrain)
    BetaML.model_save(joinpath(basefolder,"ms_xfixedpx.jld2");ms)
    xtrains  = BetaML.predict(ms,xtrain)

    ma    = BetaML.AutoEncoder(encoded_size=ae_encoded_size,layers_size=ae_hidden_layer_size,epochs=ae_base_nepochs,verbosity=BetaML.STD, cache=false)

    last_mse_val = Inf

    for n in 1:ae_max_ntrains
        verbosity >= GenFSM.STD && @info("  -- training passage: $n")
        previous_model = deepcopy(ma)
        BetaML.fit!(ma,xtrains)
        xtrainr  = BetaML.predict(ma,xtrains)
        x̂trains  = BetaML.inverse_predict(ma,xtrainr)
        x̂train   = BetaML.inverse_predict(ms,x̂trains)
        mses_train  = [BetaML.mse(xtrain[:,i],x̂train[:,i]) for i in axes(xtrain,2)]
        #means_train = mean(xtrain,dims=1)'
        xvals    = BetaML.predict(ms,xval)
        xvalsr   = BetaML.predict(ma,xvals)
        x̂vals    = BetaML.inverse_predict(ma,xvalsr)
        x̂val     = BetaML.inverse_predict(ms,x̂vals)
        mses_val = [BetaML.mse(xval[:,i],x̂val[:,i]) for i in axes(xval,2)]
        #means_val = mean(xval,dims=1)'
        mse_train_overall = StatsBase.mean(mses_train)
        mse_val_overall   = StatsBase.mean(mses_val)
        verbosity >= GenFSM.STD && @info "  -- $(n): mse_train: $(mse_train_overall)"
        verbosity >= GenFSM.STD && @info "  -- $(n): mse_val: $(mse_val_overall)"
        if mse_val_overall > last_mse_val
            verbosity >= GenFSM.STD && @info "- DONE ae pxfixed model training."
            verbosity > GenFSM.STD && @info  "  -- detailed mse_train: $(mses_train')"
            verbosity > GenFSM.STD && @info  "  -- detailed mse_val: $(mses_val')"
            ma = previous_model
            #BetaML.model_save(joinpath(basefolder,"ae_clim_m.jld2");ma)
            break
        else
            verbosity >= GenFSM.STD && @info "Further training needed..."
            last_mse_val = mse_val_overall
            if n == ae_max_ntrains
                verbosity >= GenFSM.LOW && @warn "Maximum number of training iterations reached, validation still descending. Saving the AE model, but consider further training."
                #BetaML.model_save(joinpath(basefolder,"ae_clim_m.jld2");ma)
            end
        end
    end

    xtot = BetaML.predict(ms,Matrix(xfixedpx[:,3:end]))
    xtot_reduced = BetaML.predict(ma,xtot)
    xfixedpx_reduced= hcat(xfixedpx[:,1:2],DataFrames.DataFrame(xtot_reduced,:auto))
    CSV.write(joinpath(basefolder,"xfixedpx_reduced.csv.gz"),xfixedpx_reduced;compress=true)
    BetaML.model_save(joinpath(basefolder,"ma_xfixedpx.jld2");ma)

    # Saving errors..
    x̂tots = BetaML.inverse_predict(ma,xtot_reduced)
    x̂tot  = BetaML.inverse_predict(ms,x̂tots)
    x̂fixedpx = hcat(xfixedpx[:,1:2],DataFrames.DataFrame(x̂tot,names(xfixedpx)[3:end]))

    errs = vcat([ permutedims([
        n,
        StatsBase.mean(xfixedpx[:,n]),
        BetaML.relative_mean_error(xfixedpx[:,n],x̂fixedpx[:,n]),
        BetaML.mase(xfixedpx[:,n],x̂fixedpx[:,n])
    ]) for n in names(xfixedpx)[3:end]]...,)
    errsdf = DataFrames.DataFrame(errs,["var","mean","rme","mase"])
    CSV.write(joinpath(basefolder,"xfixerpx_ae_errors.csv"),errsdf)

    return xfixedpx_reduced

end



function train_growth_model(settings, ign_growth, xclimh_reduced, xfixedpx_reduced)
    # This function trains the growth model using the IGN growth data and the autoencoded climatic and soil data
    # It returns the trained growth model.

    # creating the datasets of the growth model, for each ign data point:
    # x: forest_plot_type_dummies,xcoord,ycoord,clim_data(6),co2_concentration,pxfixed_data(6),v
    # y: dv
    # notes:
    # - dv is computed, for the observed years, as finding first i from the v2=v1(1+i)^5 eq and then multiplying it to v3
    # - the climate var are those of y2
    # - the co2 is the average of the y1:y2 years

    # -----------------------------------------------------------------------------
    # Getting options....

    verbosity = settings["verbosity"]
    verbosity >= GenFSM.STD && @info("- training growth model..")

    force_ml_train = settings["res"]["fr"]["force_ml_train"]
    basefolder     = joinpath(settings["res"]["fr"]["cache_path"],"growth_model")
    isdir(basefolder) || mkpath(basefolder)
    adj_coeff      = settings["res"]["fr"]["growth_model_training"]["vol_growth_computation_adj_coeff"]
    max_train_attempts = settings["res"]["fr"]["growth_model_training"]["max_train_attempts"]
    acceptable_rme_val = settings["res"]["fr"]["growth_model_training"]["acceptable_rme_val"]
    n_models           = settings["res"]["fr"]["growth_model_training"]["n_models"]

    # if all needed files exists and no "xclim" in force_other, then returning the saved data
    if (! ("growth" in force_ml_train)) && isfile(joinpath(basefolder,"growth_models.jld2")) && isfile(joinpath(basefolder,"growth_smodel.jld2"))
        verbosity >= GenFSM.HIGH && @info(" -- loading growth models from saved file")
        #global vol_growth_computation_parameters = JLD2.load(joinpath(basefolder,"vol_growth_computation_parameters.jld2"),"vol_growth_computation_parameters")
        return BetaML.model_load(joinpath(basefolder,"growth_smodel.jld2")), BetaML.model_load(joinpath(basefolder,"growth_models.jld2")) 
    end

    ftypes    = settings["res"]["fr"]["ftypes"]
    ohm       = BetaML.OneHotEncoder(categories = ftypes)
    co2_conc_h = CSV.read(joinpath(settings["res"]["fr"]["cache_path"],"co2_historical","co2_conc.csv"), DataFrames.DataFrame)

    nrecords = size(ign_growth,1)
    nftypes = length(ftypes)
    ae_encoded_size_cl = settings["res"]["fr"]["data_sources"]["clim"]["ae_encoded_size"]
    ae_encoded_size_px = settings["res"]["fr"]["px_fixed_data"]["ae_encoded_size"]

    # ----------------------------------------------------------------------------
    # Building the growth database (features)...

    growth_df = DataFrames.DataFrame(
        "xid" => Array{Int64,1}(undef,nrecords),
        "yid" => Array{Int64,1}(undef,nrecords),
        ["$f" => Array{Bool,1}(undef,nrecords) for f in ftypes]...,
        "x" => Array{Float64,1}(undef,nrecords),
        "y" => Array{Float64,1}(undef,nrecords),
        ["cl$(i)" => Array{Float64,1}(undef,nrecords) for i in 1:ae_encoded_size_cl ]...,
        "co2_con" => Array{Float64,1}(undef,nrecords),
        ["px$(i)" => Array{Float64,1}(undef,nrecords) for i in 1:ae_encoded_size_px ]...,
        "v" => Array{Float64,1}(undef,nrecords),
        "dv" => Array{Float64,1}(undef,nrecords),
    )


    # filling the growthdf with the data
    for (ir,r) in enumerate(eachrow(ign_growth))
        xidx = r.xid                           
        yidx = r.yid
        x = r.x
        y = r.y
        ft = vec(BetaML.predict(ohm, r.ftype))
        cl_data = xclimh_reduced[xclimh_reduced.C .== xidx .&& xclimh_reduced.R .== yidx .&& xclimh_reduced.Y .== r.y2, 4:end][1,:] # this already include the 5 years of the period considered
        pxfix_data = xfixedpx_reduced[xfixedpx_reduced.C .== xidx .&& xfixedpx_reduced.R .== yidx, 3:end][1,:]
        co2_conc = StatsBase.mean(co2_conc_h[co2_conc_h.years .>= r.y1 .&& co2_conc_h.years .<= r.y2, "co2_conc"]) # avg of the period
        v = (r.vHa2+r.vHaD+r.vHa1) / 2 # mean of the two observed periods, mortality included
        dt = r.y2 - r.y1
        
        # If I use this method I do correct in the early step of the logistic, but I go even more wrong in the late steps.
        # Using the 5 mean increments, I do overestimate early stage growth and underestimate late stages of the logistic
        # Filtering out too young plots helps keep this bias to a minimum.
        dv_rate =  ((r.vHa2 + r.vHaD)/ r.vHa1)    ^(1 / dt) - 1  # mortality included
        ydv = dv_rate * v
        # If I use this computation I got too crazy values for carrying capacity as the model doesn't bench
        # ydv = (r.vHa2 + r.vHaD - r.vHa1) / dt
        growth_df[ir,:] = [xidx, yidx,ft...,x,y,cl_data...,co2_conc,pxfix_data...,v,ydv]
    end
    CSV.write(joinpath(basefolder,"growth_df.csv"), growth_df)

    x =  Matrix(growth_df[:,3:end-1])
    y =  growth_df[:,end]

    # ----------------------------------------------------------------------------
    # Scaling...

    ms      = BetaML.Scaler(skip=1:nftypes) # skip=1:nftypes
    BetaML.fit!(ms, x)
    # hack needed because the insample co2 conc range in somehting like [350,400], but the future one goes to over 1000, if scaled only based on historical then crazy values are predicted
    co2conc_position = nftypes+ae_encoded_size_cl+3
    BetaML.parameters(ms).scalerpars.sfμ[co2conc_position]=0
    BetaML.parameters(ms).scalerpars.sfσ[co2conc_position]=0.0005
    vscp1    = BetaML.parameters(ms).scalerpars.sfμ[end]
    vscp2    = BetaML.parameters(ms).scalerpars.sfσ[end]
    vol_growth_computation_parameters = (vol_sc_par_mu=vscp1, vol_sc_par_sd=vscp2, adj_coeff=adj_coeff)
    nd = size(x,2)
    BetaML.model_save(joinpath(basefolder,"growth_smodel.jld2");ms)

    # ----------------------------------------------------------------------------
    # Model definition...

    # Nonlinear Model: NN with constraints
    # Note Logistic model:
    # - in the t/v space:   v = b/(1+c*exp(-a*t))
    # - in the dv/v space: dv = a*v -(a/b)*v^2
    #adj_coeff=0.001 
    #dvcomp(x,vscp1=vscp1,vscp2=vscp2,adj_coeff=0.001) = [x[1] * (x[3]/vscp2 - vscp1) - (adj_coeff*x[2]) * (x[3]/vscp2 - vscp1)^2]
    # JDB observation: the two parameters a,b are correlated so everytime we have higher growth rate a we have also shorted estimated carrying capacity b 
    # We should instead parametrize the logistic using the parameter Vdotmax (itself equal to ab/4) instead of a in order to get the two parameters uncorrelated. 
    # I can check with syntetic data, as in reality the network doesn't try to search for the parameters, it try to minimize the distance with observed dotV

    l1_ab   = BetaML.DenseLayer(nd-1,Int(round(1.5*(nd-1))),f=BetaML.relu)
    l1_v    = BetaML.ReplicatorLayer(1)
    l1      = BetaML.GroupedLayer([l1_ab,l1_v])
    l2_ab   = BetaML.DenseLayer(Int(round(1.5*(nd-1))),Int(round(1.5*(nd-1))),f=BetaML.relu)
    l2_v    = BetaML.ReplicatorLayer(1)
    l2      = BetaML.GroupedLayer([l2_ab,l2_v])
    l3_ab   = BetaML.DenseLayer(Int(round(1.5*(nd-1))),2,f=BetaML.relu)
    l3_v    = BetaML.ReplicatorLayer(1)
    l3      = BetaML.GroupedLayer([l3_ab,l3_v])
    l4      = BetaML.VectorFunctionLayer(3,f=Base.Fix2(GenFSM.Res_fr.vol_growth_computation,vol_growth_computation_parameters))
    layers  = [l1,l2,l3,l4]
    mgr_tpl = BetaML.NeuralNetworkEstimator(epochs=20, batch_size=16, layers=layers, fail_attempts=30, fail_epoch_ref=1, fail_epoch=4)

    partfile   = joinpath(basefolder,"growth_models_partitioned_px.csv")
    errorsfile = joinpath(basefolder,"growth_models_errors.csv")
    write(partfile, "mid,pxid,set\n")
    write(errorsfile, "mid,success,nattempts,rme_train,rme_val,rme_test,mase_train,mase_val,mase_test\n")

    # ----------------------------------------------------------------------------
    # Train function definition...

    # Train the model.. several attempts are made as sometimes (..often..) the model either doesn't converge or is linear
    # The model is considered linear if the b parameter is too high (i.e. > 2e6)
    # The model is considered converged if the relative mean error on the test set is below acceptable_rme_val
    # The model is reset and reinitialized if it doesn't converge or is linear
    # Note that the training lukely is quite bipolar: either it converge to a linear model or it converge to a "rreasonable" nonlinear model.
    function fit_hard!(mgr,xtrains,ytrain,xvals,yval;vol_growth_computation_parameters,verbosity=GenFSM.STD,max_train_attempts=30,acceptable_rme_val=0.5)
        success = false
        attempt = 1
        rme_train = Inf
        rme_val  = Inf
        # First training...
        while !success
            verbosity >= GenFSM.LOW && @info "*** Training constrained model attempt $(attempt)..."
            ŷtrain     = BetaML.fit!(mgr, xtrains, ytrain)
            ŷval       = BetaML.predict(mgr,xvals)
            rme_train  = BetaML.relative_mean_error(ytrain,ŷtrain)
            rme_val    = BetaML.relative_mean_error(yval,ŷval)
            js = rand(1:size(xtrains,1),5)
            if rme_val < acceptable_rme_val
                verbosity >= GenFSM.STD && @info "Growth model found a solution with rme_train=$(rme_train) and rme_val=$(rme_val). Let's check it is not linear..."
                linear = false
                for j in js
                    r = xtrains[j,:]
                    (aj,bj,vj,dvj) = GenFSM.Res_fr.vol_growth_computation_get_coefficients(mgr,r;adj_coeff=vol_growth_computation_parameters.adj_coeff)
                    verbosity >= GenFSM.STD && @info "b parameter: $bj"
                    if abs(bj) > 2e6 
                    linear = true
                    end
                end
                if !linear # Chacking also that the model output is not linear
                    verbosity >= GenFSM.STD && @info "Growth model is not linear, we got what we wanted !!"
                    success = true
                end
            else
                verbosity >= GenFSM.STD && @info "Growth model rme_val=$(rme_val) is not acceptable (max acceptable is $(acceptable_rme_val)). Let's try again..."
            end
            if !success 
                verbosity >= GenFSM.STD && @info "going to reset the model and try again..."
                BetaML.reset!(mgr)
                BetaML.Nn.random_init!.(BetaML.hyperparameters(mgr).layers)
            end
            if attempt > max_train_attempts
                verbosity >= GenFSM.LOW && @error "Growth model failed to find a solution after $attempt attempts."
                return (success,attempt,rme_train, rme_val)
            end
            attempt += 1
        end
        return (success,attempt,rme_train, rme_val)
    end

    # ----------------------------------------------------------------------------
    # Starting the loop over many different models...

    n_successes = 0
    for jm in 1:n_models
        #jm = 3
        verbosity >= GenFSM.LOW && @info "*********************\n*** Training j model $(jm)..."
    
        # ------------------------------------------------------------------------
        # [InLoop] Data partitioning and scaling...

        ((xtrain,xval,xtest),(ytrain,yval,ytest),(ids_train,ids_val,ids_test)) = BetaML.partition([x,y,1:size(x,1)],[0.75,0.2,0.05])
        (ntrain,nval,ntest) =  size.([xtrain,xval,xtest],1)
        CSV.write(partfile,Tables.table(hcat(fill(jm,ntrain),ids_train,fill("train",ntrain)));append=true)
        CSV.write(partfile,Tables.table(hcat(fill(jm,nval),ids_val,fill("val",nval)));append=true)
        CSV.write(partfile,Tables.table(hcat(fill(jm,ntest),ids_test,fill("test",ntest)));append=true)
        xtrains  = BetaML.predict(ms,xtrain)
        xvals    = BetaML.predict(ms,xval)
        xtests   = BetaML.predict(ms,xtest)

        # ------------------------------------------------------------------------
        # [InLoop] Model j instantialization, reset and train...
        mgr = deepcopy(mgr_tpl)
        BetaML.reset!(mgr)
        BetaML.Nn.random_init!.(BetaML.hyperparameters(mgr).layers)

        (success_flag, nattempts, rme_train, rme_val) = fit_hard!(mgr,xtrains,ytrain,xvals,yval;vol_growth_computation_parameters=vol_growth_computation_parameters,max_train_attempts=max_train_attempts,acceptable_rme_val=acceptable_rme_val) # rme_train=0.36419346513945666, rme_val=0.3733832308165264

        #rme_train=0.364483301745215, rme_val=0.37536757162340695 and rme_test=0.39376201812761435

        # ------------------------------------------------------------------------
        # [InLoop] Model error assessment & saving...
        if success_flag
            n_successes += 1
            ŷtrain     = BetaML.predict(mgr,xtrains)
            ŷval       = BetaML.predict(mgr,xvals)
            ŷtest      = BetaML.predict(mgr,xtests)
            rme_test   = BetaML.relative_mean_error(ytest,ŷtest)
            mase_train = BetaML.mase(ytrain,ŷtrain)
            mase_val   = BetaML.mase(yval,ŷval)
            mase_test  = BetaML.mase(ytest,ŷtest)
            verbosity >= GenFSM.STD && @info "Growth model $(jm) trained with rme_train=$(rme_train), rme_val=$(rme_val) and rme_test=$(rme_test). Saving the model."
            open(errorsfile, "a") do f
                write(f, "$(jm),$(success_flag),$(nattempts),$(rme_train),$(rme_val),$(rme_test),$(mase_train),$(mase_val),$(mase_test)\n")
            end
            modname = "mgr_$jm"
            overrite_model_file = (n_successes == 1) ? true : false
            BetaML.model_save(joinpath(basefolder,"growth_models.jld2"),overrite_model_file;Symbol(modname) => mgr)

        else
            @error "Growth model $(jm) training failed after $(max_train_attempts) attempts."
            (rme_train,rme_val,rme_test,mase_train,mase_val,mase_test) = (missing,missing,missing,missing,missing,missing)
        end

    end

    return BetaML.model_load(joinpath(basefolder,"growth_smodel.jld2")),BetaML.model_load(joinpath(basefolder,"growth_models.jld2")) 
end

function define_state(settings, mask)
    # This function defines the state of the region based on the prepared data
    # It is called after the growth model has been trained.
end
