module CausalTools

include("helpers.jl")
include("iptw.jl")
include("psm.jl")
include("plot.jl")

using .Helpers
using .IPTW
using .PSM
using .Plotting

const iptw = IPTW.iptw
const IPTWResult = IPTW.IPTWResult
const psm = PSM.psm
const PSMResult = PSM.PSMResult

export iptw, IPTWResult, psm, PSMResult

end
