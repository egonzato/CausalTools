using GLM, DataFrames, Statistics, StatsBase, StatsModels

function iptw(data::DataFrame,formula::FormulaTerm ;truncate::Union{Nothing, AbstractVector{Int}}=nothing,type::String="Unstabilized")
    # make sure input are in the correct format
    # build formula
    ## copy dataset
    df=copy(data)
    ## paste treatment and confounders together
    treatment=formula.lhs.sym
    # treated units
    ptrt=mean(df[!,treatment])
    # calculate probabilities
    trt_model=glm(formula,df,Binomial(),LogitLink())
    predicted=predict(trt_model)
    # define weithgs based on "type" argument
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
    return (dataset=df, model=trt_model)
end