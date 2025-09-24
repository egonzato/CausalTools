# CausalTools

CausalTools is a Julia package for causal inference, including:

* iptw: Inverse Probability of Treatment Weighting

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/egonzato/CausalTools")
```

# Load data

```
using CausalTools, CSV, DataFrames
path = joinpath(dirname(pathof(CausalTools)), "..", "data", "df.csv")
data = CSV.read(path, DataFrame)
```

# Use functions

## iptw

```
object_iptw=iptw(dataset,@formula(treatment~age+smoke+sex+cholesterol),
		 truncate=[1,99],type="Stabilized")
```

Plot weights in the two treatment groups

```
plot(object_iptw)
```
