recursive_replace!(obj,replacement_pattern) = obj
recursive_replace!(obj::AbstractArray,replacement_pattern)  = recursive_replace!.(obj::AbstractArray,Ref(replacement_pattern))
recursive_replace!(obj::AbstractString,replacement_pattern) = replace(obj,replacement_pattern)
function recursive_replace!(obj::AbstractDict,replacement_pattern)
    for (k,v) in obj
        obj[k] = recursive_replace!(v,replacement_pattern)
    end
    return obj
end

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

