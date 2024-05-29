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