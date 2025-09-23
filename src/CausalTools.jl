module CausalTools

include("helpers.jl")
include("iptw.jl")
include("plot.jl")

using .Helpers
using .IPTW
using .Plotting

const iptw = IPTW.iptw
const IPTWResult = IPTW.IPTWResult

export iptw, IPTWResult

end
