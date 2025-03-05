import GenFSM.SLOAD
import GenFSM.RES

# -----------------------------------------------------------------------------
# Testing recursive_replace!
d = Dict("a" => Dict(1 => "aaa", 2 => "bgggb"), 3=>"ccc")
SLOAD.recursive_replace!(d,"ggg" => "GGG")
@test d == Dict("a" => Dict(1 => "aaa", 2 => "bGGGb"), 3=>"ccc")

# -----------------------------------------------------------------------------
# Testing expand_dict
d_input = Dict(("a","b","c") => 10, "e" => 30, ("a","b","d") => 20)
d_output = SLOAD.expand_dict(d_input)
@test d_output == Dict("a" => Dict("b" => Dict("c" => 10, "d" => 20)), "e" => 30)
d_input  = Dict(("b",3,1) => 10, 4 => 40, ("b",3)=>30)
try
    SLOAD.expand_dict(d_input)
catch e
    @test typeof(e) == ErrorException
end

# -----------------------------------------------------------------------------
# Testing override_nested_dict!
d1  = Dict("a" => 1, "b" => Dict(1 => 10, 2 => 20, 3 => Dict(1 => 100, 2 => 200)))
d2  = Dict("b" => Dict(3 => Dict(1 => 1000), 4 => 40))
# In strict mode, the absence of key "4" in d1["b"] triggers an error:
try
    SLOAD.override_nested_dict!(d1, d2; strict=true)
catch e
    @test typeof(e) == ErrorException
end
# Reset d1 for a non-strict update:
d1  = Dict("a" => 1, "b" => Dict(1 => 10, 2 => 20, 3 => Dict(1 => 100, 2 => 200)))

# With strict=false, new keys are added:
result = SLOAD.override_nested_dict!(d1, d2; strict=false)
d1 == Dict("a" => 1, "b" => Dict(1 => 10, 2 => 20, 3 => Dict(1 => 1000, 2 => 200), 4 => 40))

## Override is a dictionray of tuples...
nested_dict = Dict("a" => Dict("b" => 1, "c" => 2), "d" => 3, "e" => Dict("f" => 4, "g" => Dict("h" => 5, "i" => 6)))
override = Dict(("e","g","h") => 50, ("e","f") => 40)
SLOAD.override_nested_dict!(nested_dict,override)
@test nested_dict["e"]["g"]["h"] == 50

## Override is a nested dictionary itself...
nested_dict = Dict("a" => Dict("b" => 1, "c" => 2), "d" => 3, "e" => Dict("f" => 4, "g" => Dict("h" => 5, "i" => 6)))
override = Dict("e" => Dict("g"=>Dict("h" => 50), "f" => 40))
SLOAD.override_nested_dict!(nested_dict,override)
@test nested_dict["e"]["g"]["h"] == 50

# -----------------------------------------------------------------------------
# Testing load_settings...
override = Dict("caller_path"=>joinpath(pwd(),"foo"),("res","fr","force_download")=>["clc"])
settings = SLOAD.load_general_settings("default","default";override=override)


@test (startswith(settings["temp_path"], "/") || settings["temp_path"][2] == ':')
@test endswith(settings["temp_path"], "/test/foo/default/temp/default")
RES.load_settings!(settings)

# Override all the other settings
SLOAD.override_nested_dict!(settings,override)

@test settings["res"]["fr"]["force_download"] == ["clc"]

default_settings = SLOAD.load_full_settings("default","default")
@test  endswith(default_settings["default_scenario_path"],"/repository/default/scenarios/default")
@test  endswith(default_settings["scenario_path"],"/repository/default/scenarios/default")
@test  endswith(default_settings["temp_path"],"/test/default/temp/default")
@test  endswith(default_settings["cache_path"],"/test/default/cache")
@test  endswith(default_settings["output_path"],"/test/default/out/default")
@test default_settings["res"]["fr"]["data_sources"]["clim"]["table_id"] == "Amon" # this shoudn't change
@test default_settings["res"]["fr"]["data_sources"]["clim"]["fixed_climate"] == true


scen_settings = SLOAD.load_full_settings("default","cc_low")
@test  endswith(scen_settings["default_scenario_path"],"/repository/default/scenarios/default")
@test  endswith(scen_settings["scenario_path"],"/repository/default/scenarios/cc_low")
@test  endswith(scen_settings["temp_path"],"/test/default/temp/cc_low")
@test  endswith(scen_settings["cache_path"],"/test/default/cache")
@test  endswith(scen_settings["output_path"],"/test/default/out/cc_low")
@test scen_settings["res"]["fr"]["data_sources"]["clim"]["fixed_climate"] == false
@test scen_settings["res"]["fr"]["data_sources"]["clim"]["table_id"] == "Amon"
@test scen_settings["res"]["fr"]["data_sources"]["clim"]["experiment_id"] == "ssp126"

overridden_settings = SLOAD.load_full_settings("default","cc_strong",override=Dict(("res","fr","data_sources","clim","table_id") =>"test", "temp_path" => "/tmp", "simulation_region" => Dict("x_lb" => 0.0) ))
@test endswith(overridden_settings["default_scenario_path"],"/repository/default/scenarios/default")
@test endswith(overridden_settings["scenario_path"],"/repository/default/scenarios/cc_strong")
@test endswith(overridden_settings["temp_path"],"/tmp")
@test endswith(overridden_settings["cache_path"],"/test/default/cache")
@test endswith(overridden_settings["output_path"],"/test/default/out/cc_strong")
@test overridden_settings["res"]["fr"]["data_sources"]["clim"]["fixed_climate"] == false
@test overridden_settings["res"]["fr"]["data_sources"]["clim"]["table_id"] == "test"
@test overridden_settings["res"]["fr"]["data_sources"]["clim"]["experiment_id"] == "ssp585"
@test overridden_settings["res"]["fr"]["data_sources"]["clim"]["member_id"] == "r1i1p1f1"
@test overridden_settings["simulation_region"]["nx"] == Int(ceil((overridden_settings["simulation_region"]["x_ub"] - overridden_settings["simulation_region"]["x_lb"]) / overridden_settings["simulation_region"]["xres"]))

