module Plotting

using Plots, StatsPlots
using ..IPTW: IPTWResult

function Plots.plot(result::IPTWResult)
    df=result.dataset
    treat_col=result.treatment
    @df df density(:weight,
                   seriestype = :density,
                   group=df[!,treat_col],
                   legend=:topright,
                   xlabel="Weight",
                   ylabel="Density",
                   title="Distribution of weights by treatment",
                   lw = 2,
                   fillrange = 0,
                   fillalpha = 0.35,)
end

end
