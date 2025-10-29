"""

"""

module psm

using DataFrames, GLM, StatsModels, StatsBase
using Statistics
using ..Helpers: check_link
using GLM: LogitLink, ProbitLink, CLogLogLink

export psm, psmResult

struct psmResult
    dataset::DataFrame
    model::StatisticalModel
    # distance, logistic / mahalabonis
    # caliper,
    # ratio
    # replacement
    # scale, on the logit or probability
    # method, greedy or optimal
end

function psm(data::DataFrame, formula::FormulaTerm; distance::Nothing="Logistic",method::Nothing="Greedy", caliper::Union{Nothing,Integer,Float64}=nothing,ratio::Integer=1,replacement::Bool=false)
    # copy data
    df=copy(data)
    # extract treatment from left side function
    treatment=formula.lhs.sym
    # make sure treatment is binary
    treatment_col=df[!,treatment]
    @assert !any(ismissing, treatment_col) "Treatment column contains missing values."
    @assert all(v -> v in (0,1), treatment_col) "Treatment must be coded as 0/1."
    # check the string given in input for type of weights

end