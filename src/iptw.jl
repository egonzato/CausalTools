module IPTW

using DataFrames, GLM, StatsModels, StatsBase
using Statistics
using ..Helpers: check_link
using GLM: LogitLink, ProbitLink, CLogLogLink

export iptw, IPTWResult

struct IPTWResult
    dataset::DataFrame
    model::StatisticalModel
    type::String
    truncate::Union{Nothing, Vector{Int}}
    treatment::Symbol
end

function iptw(data::DataFrame,formula::FormulaTerm ;truncate::Union{Nothing, AbstractVector{Int}}=nothing,type::String="Unstabilized",link::Link=LogitLink())
    # build formula
    ## copy dataset
    df=copy(data)
    ## extract treatment from left side function
    treatment=formula.lhs.sym
    # make sure treatment is binary
    treatment_col=df[!,treatment]
    vals=unique(df[!,treatment]); @assert all(v->v==0||v==1,vals) "Treatment must be coded as 0/1."
    @assert !any(ismissing,treatment_col) "Treatment column contains missing values."
    if eltype(treatment_col)<:Bool df[!,treatment]=Int.(treatment_col) end
    if eltype(treatment_col)<:AbstractFloat vals=unique(treatment_col); @assert all(v->v in (0.0,1.0),vals) "Treatment must be 0/1."; df[!,treatment]=Int.(treatment_col) end
   # check that the values given for percentiles are correct
    if truncate !== nothing
        @assert length(truncate) == 2 "Truncate must be a 2-element vector"
        @assert all(0 .<= truncate .<= 100) "Truncate must be in percentage points (0–100)"
    end
    # treated units
    ptrt=mean(df[!,treatment])
    # check that the link function given in input is among the ones available
    check_link(link)
    # calculate probabilities
    trt_model=glm(formula,df,Binomial(),link)
    predicted=predict(trt_model)
    # define weights based on "type" argument
    if type=="Unstabilized"
        weight = (df[!, treatment] ./ predicted) .+ ((1 .- df[!, treatment]) ./ (1 .- predicted))
    elseif type=="Stabilized"
        weight = ifelse.(df[!, treatment] .== 1, ptrt ./ predicted, (1 - ptrt) ./ (1 .- predicted))
    end
    # attach vectors
    ## predicted probabilities
    df.predicted=predicted
    ## weights
    df.weight=weight
    ## truncated weights
    if truncate !== nothing
        low, high = truncate
        df.weight = clamp.(df.weight, percentile(df.weight,low), percentile(df.weight,high))
    end
    # return object
    return IPTWResult(df, trt_model, type, truncate,treatment)
end

end