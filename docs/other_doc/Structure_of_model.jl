# Scratchpad

module GenFSM

export run_model

module Res

function init!!()
    res_regions = ["fr","uk"]
    println("Res init..")
    for reg in res_regions   
        init!!(Val(Symbol(reg)),1)
    end
end

end # end module Res

module Mkt

function init!!()
    mkt_regions = ["fr","de"]
    println("Mkt init..")
    for reg in mkt_regions   
        init!!(Val(Symbol(reg)),1)
    end
end

end # end module mkt

module Res_fr
import ..GenFSM.Res

function _init!!(x)
    println("Res fr")
end
function Res.init!!(::Val{:fr},x)
    Res_fr._init!!(x)
end

end

module Res_uk
import ..GenFSM.Res

function _init!!(x)
    println("Res uk")
end
function Res.init!!(::Val{:uk},x)
    Res_uk._init!!(x)
end

end

module Mkt_fr
import ..GenFSM.Mkt

function _init!!(x)
    println("mkt fr")
end
function Mkt.init!!(::Val{:fr},x)
    Mkt_fr._init!!(x)
end

end

module Mkt_de
import ..GenFSM.Mkt

function _init!!(x)
    println("mkt de")
end
function Mkt.init!!(::Val{:de},x)
    Mkt_de._init!!(x)
end

end # end Mkt_de

function run_model()
    Res.init!!()
    Mkt.init!!()
end

end # end master module GenFSM


GenFSM.run_model()

a = 1







regions = ["aa","bb"]


function foo(::Val{:aa},x)
    println("aa")
    return x
end

function foo(::Val{:bb},x)
    println("bb")
    return x
end

for reg in regions
    foo(Val(Symbol(reg)),1)
end

foo(Val(:aa),1)
foo(Val(:bb),1)