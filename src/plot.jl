module Plotting

using Plots, StatsPlots
using ..IPTW: IPTWResult

function Plots.plot(result::IPTWResult)
    df = result.dataset
    treat_col = result.treatment
    @df df density(:weight,
                     group=df[!,treat_col],
                     legend=:topright,
                     xlabel="IPTW weight",
                     ylabel="Density",
                     title="Distribution of IPTW Weights by Treatment")
end

end
