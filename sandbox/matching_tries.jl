# libraries
using Revise
using CausalTools, DataFrames, CSV, StatsModels, StatsBase
# load dataset
dataset=CSV.read("C:/Users/P095206/.julia/dev/CausalTools/data/df.CSV", DataFrame)
# use propensity score matching to create matched dataset
matcheddf=psm(dataset, @formula(treatment ~ age + smoke + cholesterol),ratio=2,distance="mahalanobis",replacement=false,caliper=Inf)
# see content
mdf=matcheddf.dataset
# number of treated present
countmap(mdf.treatment)