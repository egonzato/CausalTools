""" 
iptw(data::DataFrame,formula::FormulaTerm ;truncate::Union{Nothing, AbstractVector{Int}}=nothing,type::String="Unstabilized",link::Link=LogitLink())

Calculates inverse probability of treatment weights using logistic regression, with one of the link functions available in the GLM package. Allows to calculate stabilized and unstabilized weights, as well as applying truncation to a certain percentile of the distribution of weights.

# Arguments
    * data`::DataFrame`
    * formula`::FormulaTerm`: Treatment model to be used in the GLM regression, use @formula() to create it.
    * truncate`::Union{Nothing, AbstractVector{Int}}`: Indicates whether truncation to weights of the untreated units should be applied. Indicate here the percentiles at which the distribution should be truncated. Default is set to "nothing".
    * type`::String="Unstabilized"`: String that defines whether stabilized or unstabilized weights should be calculated. Possible strings for this arguments are "Stabilized" and "Unstabilized". The latter is the default.
    * link`::Link=LogitLink()`: Link function to be used in GLM. Available functions are `LogitLink()`, which is the default, or `ProbitLink()`.

# Output
Returns a `IPTWResult` struct, which contains:
    * `dataset::DataFrame`: Original dataframe with two new columns, one with treatment probability and one with the weights calculated.
    * `model::StatisticalModel`: GLM model computed with the options specified.
    * `type::String`: String with type of weights requested.
    * `truncate::Union{Nothing, Vector{Int}}`: Percentiles used to perform truncation, if defined.
    * `treatment::Symbol`: Treatment vector.

# Example

```julia-repl
julia> object_iptw=iptw(dataset,@formula(treatment~age+smoke+sex+cholesterol),truncate=[1,99],type="Stabilized")

CausalTools.IPTW.IPTWResult(1000×8 DataFrame)
  Row │ smoke  sex    age      cholesterol  treatment  blood_pressure  predicted  weight   
      │ Bool   Bool   Float64  Float64      Bool       Float64         Float64    Float64
──────┼────────────────────────────────────────────────────────────────────────────────────
    1 │ false  false  46.7905      203.97       false         140.774   0.242294  0.848614
    2 │ false   true  62.1994      194.793       true         133.088   0.271616  1.31436
    3 │  true  false  38.7231      179.85       false         120.273   0.444942  1.15844
  ⋮   │   ⋮      ⋮       ⋮          ⋮           ⋮            ⋮             ⋮         ⋮
  999 │ false   true  43.8346      228.986      false         116.383   0.362175  1.00811
 1000 │ false   true  59.0516      191.917      false         130.752   0.281389  0.894782


plot(object_iptw)

wdataset=object_iptw.dataset

```

"""
module IPTW

using DataFrames, GLM, StatsModels, StatsBase
using Statistics
using ..Helpers: check_link
using GLM

export iptw, IPTWResult

struct IPTWResult
    dataset::DataFrame
    model::StatisticalModel
    type::String
    truncate::Union{Nothing, Vector{Int}}
    treatment::Symbol
end

function iptw(data::DataFrame,formula::FormulaTerm; truncate::Union{Nothing, AbstractVector{Int}}=nothing,type::String="Unstabilized",link::Link=LogitLink())
    ## copy dataset
    df=copy(data)
    ## extract treatment from left side function
    treatment=formula.lhs.sym
    # make sure treatment is binary
    treatment_col=df[!,treatment]
    @assert !any(ismissing, treatment_col) "Treatment column contains missing values."
    @assert all(v -> v in (0,1), treatment_col) "Treatment must be coded as 0/1."
    # check the string given in input for type of weights
    @assert type == "Stabilized" || type == "Unstabilized" "Invalid type argument: must be 'Stabilized' or 'Unstabilized'"
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