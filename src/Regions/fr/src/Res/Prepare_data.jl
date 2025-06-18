# Prepare the data for the French Rres_fr module.

function prepare_data!(settings, mask)
    settings["verbosity"] >= LOW && @info("Preparing res data for the French region.")
    # This function prepares the data for the French region.
    # It is called after the data has been downloaded and layers saved as tiff.
    
    ign_state, ign_growth = prepare_ign_data!(settings)
    # the ae_clim and ae_soil models could be used as sort of pretrain and then chained vertically together to form the final growth/mortality model instead of using only the reduced form, but it would takes a lot of computational power in prediction, that is not suitable 
    xclimh_reduced, scaler_clim_m, ae_clim_m  = train_autoencode_clim(settings, mask, ign_growth)
    xclimf_reduced   = GenFSM.Res_fr.predict_autoencoder_clim(settings, mask, scaler_clim_m, ae_clim_m)
    scaler_clim_m, ae_clim_m = nothing, nothing # removing the climatic ae models from memory, as they are not needed anymore
    
    xsoil_reduced       = train_autoencode_soil!(settings, mask)
    #res_growth_m    = train_growth_model(settings, ign_state ign_growth, xclimh_reduced, xsoil_reduced)
    #res_mortality_m = train_mortality_model(settings, ign_state, ign_growth, xclimh_reduced, xsoil_reduced)
    #define_state(settings,mask)
    # xclimf_reduced                  = predict_autoencoder_clim(settings, mask, scaler_clim_m, ae_clim_m)
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
y: year
vHa: volume per hectare
d: avg diameter
spg: species group BR,MIX,CON
gov: government HF,MIX,COP

"""
function prepare_ign_data(settings)

    # This function prepares the IGN data for the French region
    settings["verbosity"] >= STD && @info("- preparing IGN data for growth/mortality model training and set initial forest resource condition")

    force_other = settings["res"]["fr"]["force_other"]
    force_ml_train = settings["res"]["fr"]["force_ml_train"]
    basefolder  = joinpath(settings["res"]["fr"]["cache_path"],"forinv")

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

    map_espar = Dict(species.spcode .=> species.group)
    # MAnually add codes that are not in the species list....
    map_espar["332G"] = "br"  # populus
    map_espar["29AF"] = "br"  # other broadleaves
    map_espar["68CE"] = "con" # other coniferous 
    map_espar["25E3"] = "br"  # some salix

    pids  = unique(pobs.IDP)
    years = unique(pobs.CAMPAGNE)

    map_structure_sfo   = Dict(0 => "", 1 => "HF", 2 => "", 3 => "MIX", 4 => "COP")  # skipping irregular HF (2) and no structure (0)
    map_structure_sver  = Dict("0" => "", "X" => "", "2" => "HF", "3" => "COP", "4" => "", "5" => "MIX", "6" => "HF") # skipping no structure/cutted (0,X) and irregular HF (4)
    map_recent_cuts = Dict(0 => false, 1 => true, 2 => true)

    function vHaContribution(v,c13)
        (ismissing(v) || ismissing(c13)) && return missing
        if c13 < 0.705
            return v/(6^2*pi/(100*100))
        elseif c13 < 1.175
            return v/(9^2*pi/(100*100))
        else 
            return v/(15^2*pi/(100*100))
        end
    end
    train_fields = ["ESPAR","CAMPAGNE","SER","C13","HTOT","V"]
    v_cidx   = findfirst(==("V"), train_fields)
    
    if !("tree_volumes" in force_ml_train) && (isfile(joinpath(basefolder,"tree_imputer_model.jd")))
        tree_imputer_model = BetaML.model_load(joinpath(basefolder,"tree_imputer_model.jd"),"tree_imputer_model")
    else
        # Option 1: do the imputation and save the model...
        @info("Training individual trees imputation model...")

        train_tpobs = DataFrames.dropmissing(tpobs[:,train_fields])
        (N,D) = size(train_tpobs)
        Nsample = nsample_v_imputation
        sids = StatsBase.sample(axes(train_tpobs, 1), Nsample; replace = false, ordered = false)
        training_tpobs = view(train_tpobs, sids, train_fields)
        tree_imputer_model = BetaML.RandomForestEstimator(verbosity=BetaML.FULL,cache=false, n_trees=50, oob=true)

        BetaML.fit!(tree_imputer_model, Matrix(training_tpobs[:,1:end-1]), training_tpobs[:,end])
        # saving model
        BetaML.model_save(joinpath(basefolder,"tree_imputer_model.jd");tree_imputer_model=tree_imputer_model)
        @info("...done individual trees imputation model")
        @info(BetaML.info(tree_imputer_model))
    end

    if !( "prediction_hist_tree_volumes" in force_other) && (isfile(joinpath(basefolder,"tpobs.csv")))
        tpobs = CSV.File(joinpath(basefolder,"tpobs.csv")) |> DataFrames.DataFrame
    else
        @info("Starting predicting individual trees volumes...")
        x_topredict = copy(tpobs[:,train_fields[1:end-1]])
        vest = BetaML.predict(tree_imputer_model, Matrix(x_topredict))
        tpôbs      = deepcopy(tpobs)
        tpôbs.V    = vest
        # saving data data
        CSV.write(joinpath(basefolder,"tpobs.csv"),tpôbs)
        tpobs = tpôbs
    end

    # Computing, after imputation, the vHa for each tree
    tpobs.vHa  = vHaContribution.(tpobs.V,tpobs.C13)

    # Selecting only points with 2 visits
    cmap = StatsBase.countmap(pobs.IDP)
    pids2 = filter(x -> cmap[x] == 2, pids)
    #pobs2  = pobs[in.(pobs.IDP, Ref(pids2)), :]

    # Step 1: computing dfgrowth data
    @info("Starting building dfgrowth DataFrame...")
    counts = fill(0,7)
    dfgrowth = DataFrames.DataFrame(idp = Int64[], x = Float64[], y = Float64[], y1 = Int64[], y2 = Int64[],
                    vHa1 = Float64[], vHa2 = Float64[], vHaD = Float64[], d1 = Float64[], d2 = Float64[],
                    spg = String[], gov = String[])
    for (i,p) in enumerate(pids2)
        counts[1] += 1
        i%1000 == 0 && println("p: $p np: $(length(pids2)) counts: $counts")
        #p = pids2[end]
        #p = 146283 
        mytrees2 = tpobs[tpobs.IDP .== p .&& tpobs.VISITE .==2, :]
        size(mytrees2,1) < 2        && continue  # skip points with < 2 (dead or alive) trees on the second visit
        counts[2] += 1
        mytrees2.vegstate = [ismissing(t.VEGET5) ? t.VEGET : t.VEGET5  for t in eachrow(mytrees2)] # some second visit trees have veg state in VEGET instead of VEGET5 :-()
        mytrees2_cutted = mytrees2[in.(mytrees2.vegstate, Ref(["6","7"])), :]
        size(mytrees2_cutted,1) == 0 || continue  # skip points with cuts
        counts[3] += 1
        visit1 = pobs[pobs.IDP .== p .&& pobs.VISITE .== 1, :][1,:]
        # Vertical structure
        vs= ""
        if !ismissing(visit1.SFO)
            vs = map_structure_sfo[visit1.SFO]
        elseif !ismissing(visit1.SVER) 
            vs = map_structure_sver[visit1.SVER]
        end
        vs !== "" || continue # skip points with no vertical structure
        counts[4] += 1
        mytrees2_alive  = mytrees2[mytrees2.vegstate .== "0", :]
        mytrees1 = tpobs[tpobs.IDP .== p .&& tpobs.VISITE .==1, :]
        # Setting the vHa for the alive trees on the second visit as minima as their vHa on the first visit, if present
        for (j,t) in enumerate(eachrow(mytrees2_alive))
            t1s =  mytrees1[mytrees1.A .== t.A,:]
            size(t1s,1) == 0 && continue    # tree not present on first visit, it's new mortality
            t1 = t1s[1,:]
            ismissing(t1.vHa) && continue # missing volume on first visit tree
            t.vHa = max(t.vHa, t1.vHa) # setting the vHa for the alive trees on the second visit as minima as their vHa on the first visit
        end
        mytrees2_alive.vHa = vHaContribution.(mytrees2_alive.V, mytrees2_alive.C13) # computing vHa for alive trees on second visit
        vHa2 = sum(mytrees2_alive.vHa)
        ismissing(vHa2) && continue # skip points with no volume info on alive trees on second visit 
        counts[5] += 1
        # ended filters... (not really..)
        mytrees1.vegstate = [ismissing(t.VEGET5) ? t.VEGET : t.VEGET5  for t in eachrow(mytrees1)] 
        mytrees1_alive  = mytrees1[mytrees1.vegstate .== "0", :]
        mytrees1_alive.vHa = vHaContribution.(mytrees1_alive.V, mytrees1_alive.C13) # computing vHa for alive trees on first visit
        visit2 = pobs[pobs.IDP .== p .&& pobs.VISITE .== 2, :][1,:]
        mytrees2_dead      = mytrees2[.! in.(mytrees2.vegstate, Ref(["0","6","7"])), :]
        mytrees2_dead.vHA .= 0.0
        countd1 = size(mytrees2_dead,1)
        # Removing dead tress that were already dead on the first visit and assigning to deatch trees (that have no measured diameter) the vHa of the first visit
        to_keep = fill(true,size(mytrees2_dead,1))
        for (j,t) in enumerate(eachrow(mytrees2_dead))
            t1s =  mytrees1[mytrees1.A .== t.A,:]
            size(t1s,1) == 0 && continue    # tree not present on first visit, it's new mortality
            t1 = t1s[1,:]
            ismissing(t1.VEGET) && continue # not reported state on visit 1, assuming alive, so new mortality
            (t1.VEGET .== "0") || (to_keep[j] = false) # was already not alive 
            t.vHa = vHaContribution(t1.V,t1.C13) # computing vHa for dead trees on second visit
        end
        mytrees2_dead = mytrees2_dead[to_keep,:]
        countd2 = size(mytrees2_dead,1)
        # Removing points with trees that don't have hTOT, as their volumes are not reliable
        check = false 
        for t in eachrow(vcat(mytrees1_alive, mytrees2_alive))
            ismissing(t.HTOT) && (check = true; break)
        end
        # if I enable this one I have no more points :-|
        #check == false || continue # skip points with no HTOT on alive trees
        counts[6] += 1
        # Defining the species group
        mytrees1_alive.spgr = [map_espar[t.ESPAR] for t in eachrow(mytrees1_alive)]
        vHaBr  = sum( mytrees1_alive[mytrees1_alive.spgr .== "br", "vHa"])
        vHaCon = sum( mytrees1_alive[mytrees1_alive.spgr .== "con", "vHa"])
        vHaTot = vHaBr .+ vHaCon
        vHaBrRatio = vHaBr ./ vHaTot
        if vHaBrRatio < 0.1 
            spgr = "con"
        elseif vHaBrRatio < 0.9
            spgr = "mixsp"
        else
            spgr = "br"
        end

        # Removing points with no alive trees on first visit and both no new death trees and no alive trees on visit 2
        size(mytrees1_alive,1) == 0 && size(mytrees2_dead,1) == 0 && size(mytrees2_alive,1) == 0 && continue
        counts[7] += 1

        vHa1 = sum(mytrees1_alive.vHa)
        vHaD = sum(mytrees2_dead.vHa)
        d1 = (size(mytrees1_alive,1) == 0) ? 0.0 : StatsBase.mean(mytrees1_alive.C13)  # some plots have no trees on the first visit, so no diameter. This is ok when the trees are too small to be measured
        d2 = (size(mytrees2_alive,1) == 0) ? d1 : StatsBase.mean(mytrees2_alive.C13)
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
        # Now I have everything I need to create the row for the growing_points DataFrame
        push!(dfgrowth, [p, visit1.XL, visit1.YL, visit1.CAMPAGNE, visit2.CAMPAGNE, vHa1, vHa2, vHaD, d1, d2, spgr, vs])
    end
    dfgrowth.dV     = (dfgrowth.vHa2+dfgrowth.vHaD) .- dfgrowth.vHa1
    dfgrowth.dV2    = dfgrowth.vHa2 .- dfgrowth.vHa1
    dfgrowth.mShare = dfgrowth.vHaD ./ (dfgrowth.vHa2 .+ dfgrowth.vHaD) # share of dead volume in the total volume change
    CSV.write(joinpath(basefolder,"dfgrowth.csv"),dfgrowth)

    @info("Done building dfgrowth DataFrame. Counts: $counts")

    # Step 2: building dfstate DataFrame
    @info("Starting building dfstate DataFrame...")
    pid_latest = pobs[ in.(pobs.CAMPAGNE, Ref(maximum(pobs.CAMPAGNE)-2 : maximum(pobs.CAMPAGNE))), "IDP"]
    #443593 in first visit and 218916 in second visit

    counts = fill(0,4)
    dfstate = DataFrames.DataFrame(idp = Int64[], x = Float64[], y = Float64[], year=Int64[], vHa = Float64[], d = Float64[], spg = String[], gov = String[])
    for (i,p) in enumerate(pid_latest)
        counts[1] += 1
        i%1000 == 0 && println("p: $p np: $(length(pid_latest)) counts: $counts")
        #p = pid_latest[end]
        #p = 146283 
        visit = pobs[pobs.IDP .== p, :][1,:]

        in.(visit.CSA, Ref(["1","3","5","2"])) || continue # skip points with no forest structure (open, closed forest or popplars)
        counts[2] += 1

        mytrees = tpobs[tpobs.IDP .== p, :]
        mytrees.vegstate = [ismissing(t.VEGET5) ? t.VEGET : t.VEGET5  for t in eachrow(mytrees)] # some second visit trees have veg state in VEGET instead of VEGET5 :-()
        mytrees_alive  = mytrees[mytrees.vegstate .== "0", :]
        mytrees_alive.vHa = vHaContribution.(mytrees_alive.V, mytrees_alive.C13) # computing vHa for alive trees on second visit
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
        push!(dfstate, [p, visit.XL, visit.YL, visit.CAMPAGNE, vHa, d, spgr, vs])
    end
    @info("Done building dfgrowth DataFrame. Counts: $counts")

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
    nxclimh   = nC*nR*(length(hyears)-ae_nyears+1)
  
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
    nxclimf  = nC*nR*(length(fyears))
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


function train_growth_model(settings, ign_growth, ae_clim_m, ae_soil_m)
    # This function trains the growth model using the IGN growth data and the autoencoded climatic and soil data
    # It returns the trained growth model.
    return nothing
end

function define_state(settings, mask)
    # This function defines the state of the region based on the prepared data
    # It is called after the growth model has been trained.
end
