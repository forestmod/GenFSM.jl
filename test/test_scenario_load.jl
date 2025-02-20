import GenFSM.SLOAD

nested_dict = Dict("a" => Dict("b" => 1, "c" => 2), "d" => 3, "e" => Dict("f" => 4, "g" => Dict("h" => 5, "i" => 6)))
override = Dict(("e","g","h") => 50, ("e","f") => 40)
SLOAD.override_nested_dict!(nested_dict,override)

@test nested_dict["e"]["g"]["h"] == 50

settings = SLOAD.load_settings("default","default";override=Dict("caller_path"=>joinpath(pwd(),"foo"),("res","fr","force_download")=>["clc"]))

@test (startswith(settings["temp_path"], "/") || settings["temp_path"][2] == ':')
@test endswith(settings["temp_path"], "/test/foo/temp")

@test settings["res"]["fr"]["force_download"] == ["clc"]