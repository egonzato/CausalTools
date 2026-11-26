# libraries
using Revise
using Pkg, CausalTools, DataFrames, CSV, StatsModels, StatsBase
# load dataset
dataset=CSV.read(joinpath(@__DIR__, "..", "data", "df.csv"), DataFrame)
# use propensity score matching to create matched dataset
matcheddf=psm(dataset, @formula(treatment ~ age + smoke + cholesterol),ratio=2,distance="mahalanobis",replacement=false,caliper=Inf)
# see content
mdf=matcheddf.dataset
# number of treated present
countmap(mdf.treatment)