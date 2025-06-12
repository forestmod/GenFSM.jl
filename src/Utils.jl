

"""

$(TYPEDEF)

Many models and functions accept a `verbosity` parameter.

Choose between: `NONE`, `LOW`, `STD` [default], `HIGH` and `FULL`.
"""
@enum Verbosity NONE=0 LOW=10 STD=20 HIGH=30 FULL=40

verbosity_map = Dict("NONE" => NONE, "LOW" => LOW, "STD" => STD, "HIGH" => HIGH, "FULL" => FULL)