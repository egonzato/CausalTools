
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
data = CSV.read("lalonde.csv", DataFrame)
```

# Use functions

## iptw

```
iptw_obj = iptw(data, @formula(treat ~ age + race), truncate=[1,99], type="Stabilized")
```

Plot weights in the two treatment groups

```
plot(iptw_obj)
```
