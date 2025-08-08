

"""

$(TYPEDEF)

Many models and functions accept a `verbosity` parameter.

Choose between: `NONE`, `LOW`, `STD` [default], `HIGH` and `FULL`.
"""
@enum Verbosity NONE=0 LOW=10 STD=20 HIGH=30 FULL=40

verbosity_map = Dict("NONE" => NONE, "LOW" => LOW, "STD" => STD, "HIGH" => HIGH, "FULL" => FULL)




"""
    df_to_dict(df, dim_col, value_cols)
    Convert a DataFrame to a dictionary specifying the column(s) to use  as keys and the columns to use as values.
"""
function to_dict(df::DataFrames.AbstractDataFrame, dim_cols, value_cols)
    # normalize inputs
    dim_syms = isa(dim_cols, AbstractVector)   ? Symbol.(dim_cols)   : [Symbol(dim_cols)]
    val_syms = isa(value_cols, AbstractVector) ? Symbol.(value_cols) : [Symbol(value_cols)]
    to_return = Dict{Tuple, Vector}()
    for r in DataFrames.eachrow(df)
        key_labels = []
        [push!(key_labels,r[d]) for d in dim_syms]
        to_return[(key_labels...,)] = Vector(r[val_syms])
    end
    return to_return
end