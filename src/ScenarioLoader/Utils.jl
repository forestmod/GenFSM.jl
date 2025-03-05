
"""
    recursive_replace!(dict,replacement_pattern)

Replace a pattern in the values a possibly nested dictionary

```julia
d = Dict("a" => Dict(1 => "aaa", 2 => "bgggb"), 3=>"ccc")
recursive_replace!(d,"ggg" => "GGG") # Dict("a" => Dict(1 => "aaa", 2 => "bGGGb"), 3=>"ccc")
```

"""
recursive_replace!(obj,replacement_pattern) = obj
recursive_replace!(obj::AbstractArray,replacement_pattern)  = recursive_replace!.(obj::AbstractArray,Ref(replacement_pattern))
recursive_replace!(obj::AbstractString,replacement_pattern) = replace(obj,replacement_pattern)
function recursive_replace!(obj::AbstractDict,replacement_pattern)
    for (k,v) in obj
        obj[k] = recursive_replace!(v,replacement_pattern)
    end
    return obj
end

# Helper function that recursively inserts a value into the dictionary
# following the structure defined by the tuple of keys.
function insert_nested!(d::AbstractDict, keys::Tuple, value)
    if length(keys) == 1
        # Base case: if the tuple has only one element, assign the value.
        d[keys[1]] = value
    else
        # Use the first element of the tuple as the key.
        k = keys[1]
        # If the key does not exist, create a new nested dictionary.
        if !haskey(d, k)
            d[k] = Dict{Any,Any}()
        # If the key exists but is not a dictionary, report a conflict.
        elseif !isa(d[k], Dict)
            error("Key conflict: expected a dictionary at key $k")
        end
        # Recursively insert into the nested dictionary using the rest of the tuple.
        insert_nested!(d[k], keys[2:end], value)
    end
end

"""
   expand_dict(d::Dict)

Expand dictionary keys that are tuples into nested dictionaries.

# Example

```julia
d_input = Dict(("a","b","c") => 10, "e" => 30, ("a","b","d") => 20)
d_output = expand_dict(d_input) # Dict("a" => Dict("b" => Dict("c" => 10, "d" => 20)), "e" => 30)
```

"""
function expand_dict(d::AbstractDict)
    result = Dict{Any,Any}()
    for (k, v) in d
        if isa(k, Tuple)
            insert_nested!(result, k, v)
        else
            # If the key is not a tuple, simply assign it.
            result[k] = v
        end
    end
    return result
end

"""
    override_nested_dict!(d1::AbstractDict, d2::AbstractDict; strict::Bool=true, expand=true)

Override the values in the dictionary `d1` with the values in the dictionary `d2`.

d2 can be either another nested dictionary or a dictionary where the different levels are represented as tuple keys, e.g. `d2 = Dict("a"=>Dict("b"=>1))` or `d2 = Dict(("a","b")=>1)`.
If d2 has new keys and `strict` is enabled, an error is raised.

"""
function override_nested_dict!(d1::AbstractDict, d2::AbstractDict; strict::Bool=true, expand=true)
    if expand # only once
        d2 = expand_dict(d2) # if the keys of the dictionary with the values that override is expressed as tuple, we transform the dictionary in a nested dictionary
    end
    for (k, v) in d2
        if haskey(d1, k)
            # If both corresponding values are dictionaries, recurse.
            if isa(d1[k], Dict) && isa(v, Dict)
                override_nested_dict!(d1[k], v; strict=strict, expand=false)
            else
                # Otherwise, simply override the value.
                d1[k] = v
            end
        else
            if strict
                error("Key $(k) does not exist in the main dictionary.")
            else
                d1[k] = v
            end
        end
    end
    return d1
end

#=
function override_nested_dict!(nested_dict,override_elements)
    for o_element in override_elements
        td = nested_dict
        o_ks = o_element[1]; o_v = o_element[2];
        (typeof(o_ks) <: Tuple) || (o_ks = (o_ks,))
        for (ik,k) in enumerate(o_ks)
            k in keys(td) || error("Override error: Key $k not found in container $td")
            ik == length(o_ks) ?  (td[k] = o_v)  :  (td = td[k]) 
        end
    end
end

=#