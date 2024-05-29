
# ------------------------------------------------------------------------------
println("*** Testing `override_nested_dict!`...")
a_nested_dict     = Dict("a"=>1,"b"=>Dict("c"=>2,"d"=>[Dict("e"=>2),Dict("f"=>3)]),"g"=>[1,"h"])
override_elements = ["a"=>10,("g",2)=>"l",("b","d",2)=>Dict("i"=>4)]
override_nested_dict!(a_nested_dict,override_elements)
@test a_nested_dict == Dict("a"=>10,"b"=>Dict("c"=>2,"d"=>[Dict("e"=>2),Dict("i"=>4)]),"g"=>[1,"l"])

a_nested_dict     = Dict("a"=>1,"b"=>Dict("c"=>2,"d"=>[Dict("e"=>2),Dict("f"=>3)]),"g"=>[1,"h"])
override_elements = Dict("a"=>10,("g",2)=>"l",("b","d",2)=>Dict("i"=>4))
override_nested_dict!(a_nested_dict,override_elements)
@test a_nested_dict == Dict("a"=>10,"b"=>Dict("c"=>2,"d"=>[Dict("e"=>2),Dict("i"=>4)]),"g"=>[1,"l"])

a_nested_dict     = Dict("a"=>1,"b"=>Dict("c"=>2,"d"=>[Dict("e"=>2),Dict("f"=>3)]),"g"=>[1,"h"])
override_elements = []
override_nested_dict!(a_nested_dict,override_elements)
@test a_nested_dict == Dict("a"=>1,"b"=>Dict("c"=>2,"d"=>[Dict("e"=>2),Dict("f"=>3)]),"g"=>[1,"h"])
