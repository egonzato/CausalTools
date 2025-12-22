"""
psm(data::DataFrame, formula::FormulaTerm; distance::String="probability", ratio::Int=1, replacement::Bool=false, caliper::Float64=Inf, order::String="largest", seed::Int=181299, discard::Bool=false)

Performs propensity score matching (PSM) or Mahalanobis distance matching using nearest neighbor algorithm. Matches treated units to control units based on specified distance metric.

# Arguments
* `data::DataFrame`: Input dataset containing treatment and covariates.
* `formula::FormulaTerm`: Formula specifying treatment ~ covariates, use `@formula()` to create it.
* `distance::String="probability"`: Distance metric. Options: `"probability"` (raw PS), `"logit"` (logit-transformed PS), `"mahalanobis"` (Mahalanobis on covariates).
* `ratio::Int=1`: Number of controls matched per treated unit.
* `replacement::Bool=false`: Whether to allow control units to be matched multiple times.
* `caliper::Float64=Inf`: Caliper width in standard deviations of distance metric. `Inf` = no caliper.
* `order::String="largest"`: Ordering of treated units for sequential matching. Options: `"largest"`, `"smallest"`, `"random"`.
* `seed::Int=181299`: Random seed for `"random"` ordering.
* `discard::Bool=false`: Whether to discard treated units without complete matches (`true`) or keep partial matches (`false`).

# Output
Returns a `PSMResult` struct containing:
* `dataset::DataFrame`: Matched dataset with original columns + `id`, `cluster` (matching groups).
* `unmatched::Vector{Int}`: IDs of unmatched treated units.
* `distance::String`: Distance metric used.
* `ratio::Int`: Matching ratio used.
* `replacement::Bool`: Replacement setting.
* `n_matched::Int`: Number of matched observations.

# Example

"""
module PSM

using DataFrames, StatsModels, GLM, StatsFuns, Random, StatsBase, LinearAlgebra

export psm, PSMResult

struct PSMResult
    dataset::DataFrame
    unmatched::Vector{Int}
    distance::String
    ratio::Int
    replacement::Bool
end
# probabilities
function calculate_probabilities(dataset, formula)
    model = glm(formula, dataset, Binomial())
    dataset.distance = predict(model)
    return dataset
end
# logit
function calculate_logit(dataset, formula)
    model = glm(formula, dataset, Binomial())
    dataset.distance = logit.(predict(model))
    return dataset
end
# mahalanobis
function calculate_mahalanobis(x, Y, Σ_inv)
    n = size(Y, 1)
    dists = zeros(n)
    for i in 1:n
        x_row = view(Y, i, :)
        diff_i = x_row - x
        dists[i] = sqrt(diff_i' * Σ_inv * diff_i)
    end
    return dists
end


function psm(dataset::DataFrame, formula::FormulaTerm; distance="probability",ratio=1,replacement=false,caliper::Float64=Inf,order="largest",seed=181299,discard=false)
    # Extract treatment column
    treatment=formula.lhs.sym
    # warnings
      # ratio must be an integer ≥ 1
    if !(isa(ratio, Integer) && ratio ≥ 1)
        throw(ArgumentError("`ratio` must be an integer ≥ 1. Got: $ratio"))
    end
    # replacement must be true or false
    if !(replacement === true || replacement === false)
        throw(ArgumentError("`replacement` must be true or false. Got: $replacement"))
    end
    # discard must be true or false
    if !(discard === true || discard === false)
        throw(ArgumentError("`discard` must be true or false. Got: $discard"))
    end
    # distance must be "logit" or "probability"
    allowed_distances = ["logit", "probability","mahalanobis"]
    if !(distance in allowed_distances)
        throw(ArgumentError("`distance` must be one of $allowed_distances. Got: $distance"))
    end
    # caliper must be ≥ 0
    if !(isa(caliper, Real) && caliper ≥ 0)
        throw(ArgumentError("`caliper` must be ≥ 0. Got: $caliper"))
    end
    # order must be "largest", "smallest" or "random"
    allowed_orders = ["largest", "smallest", "random"]
    if !(order in allowed_orders)
        throw(ArgumentError("`order` must be \"largest\", \"smallest\", or \"random\". Got: $order"))
    end
    # seed must be an integer
    if !(isa(seed, Integer))
        throw(ArgumentError("`seed` must be an integer. Got: $seed"))
    end
    # check forid column
    if :id ∉ names(dataset)
        dataset.id = 1:nrow(dataset)
    end
    # treatment column exists
    #if treatment ∉ names(dataset)
    #    throw(ArgumentError("Treatment '$treatment' not found in dataset"))
    #end
    # covariates exist
    if distance == "mahalanobis"
        covars = [t.sym for t in formula.rhs]
        missing_covars = setdiff(covars, names(dataset))
        if !isempty(missing_covars)
            throw(ArgumentError("Mahalanobis covariates missing: $missing_covars"))
        end
    end
    # copy dataset
    df=deepcopy(dataset)
    df.id = 1:nrow(df)
    # compute distances and sort descending
    if distance=="probability"
        df=calculate_probabilities(df, formula)
    elseif distance=="logit"
        df=calculate_logit(df, formula)
    end
    # order dataset given a preference, if selected either probability or logit
    if distance in allowed_distances[1:2]
        if order=="largest"
            sort!(df,:distance,rev=true)
        elseif order=="smallest"
            sort!(df,:distance)
        elseif order=="random"
            Random.seed!(seed)
            df=df[shuffle(1:nrow(df)),:]
        end
    end
    # extract treatment column
    treatment_col = df[!, treatment]
    # check treatment column
    unique_vals = unique(treatment_col)
    if !all(x -> x === true || x === false, unique_vals)
        throw(ArgumentError("The treatment variable '$treatment' must contain only true or false values. Found: $(unique_vals)"))
    end
    # filter out treated and untreated
    untrt=df[treatment_col.==false,:]
    trt=df[treatment_col.!=false,:]
    # non-empty pools
    if nrow(trt) == 0
        throw(ArgumentError("No treated units found"))
    end
    if nrow(untrt) == 0  
        throw(ArgumentError("No control units found"))
    end
    # define std caliper
    if distance=="probability"
        std_caliper=caliper * std(logit.(df.distance))
    elseif distance=="logit"
        std_caliper=caliper * std((df.distance))
    end
    # initialize matched dataset and used controls
    matched=DataFrame(id=Int[], cluster=Int[])
    unmatched=[]
    pair=DataFrame(id=Int[], cluster=Int[])
    # search controls
    if distance in allowed_distances[1:2]
        if !replacement
            for r in 1:ratio
                for i in 1:nrow(trt)
                    # define filtered at the beginning
                    filtered=copy(untrt)
                    # stop if ran out of controls
                    if nrow(untrt)==0
                        @warn("Ran out of controls")
                        break
                    end
                    # define filtered dataset if caliper exists
                    if caliper!=Inf && distance=="probability"
                        filtered.logit=logit.(filtered.distance)
                        filtered.diff=abs.(filtered.logit .- logit.(trt.distance[i])) 
                        filtered=filtered[filtered.diff .<= std_caliper,:]
                    elseif caliper!=Inf && distance=="logit"
                        filtered.logit=filtered.distance
                        filtered.diff=abs.(filtered.logit .- trt.distance[i]) 
                        filtered=filtered[filtered.diff .<= std_caliper,:]
                    else 
                        filtered.diff=abs.(filtered.distance .- trt.distance[i]) 
                    end
                    # sort by difference to get the two closest
                    sort!(filtered,[:diff])
                    # skip if no control is available
                    if nrow(filtered)==0
                        # add to unmatched
                        push!(unmatched,trt.id[i])
                        continue
                    end
                    # append only control if iteration is not the first
                    if r == 1
                        pair=DataFrame(id=[trt.id[i];filtered.id[1]],cluster=fill(i,2))
                    else 
                        pair=DataFrame(id=[filtered.id[1]],cluster=[i])
                    end
                    # push into matched
                    append!(matched,pair)
                    # delete controls that have been used
                    untrt=untrt[.!in.(untrt.id, Ref(filtered.id[1])), :]
                end 
            end
        else
            for i in 1:nrow(trt)
                # define filtered at the beginning
                filtered=copy(untrt)
                # define filtered dataset if caliper exists
                if caliper!=Inf && distance=="probability"
                    filtered.logit=logit.(filtered.distance)
                    filtered.diff=abs.(filtered.logit .- logit.(trt.distance[i])) 
                    filtered=filtered[filtered.diff .<= std_caliper,:]
                elseif caliper!=Inf && distance=="logit"
                    filtered.logit=filtered.distance
                    filtered.diff=abs.(filtered.logit .- trt.distance[i]) 
                    filtered=filtered[filtered.diff .<= std_caliper,:]
                else 
                    filtered.diff=abs.(filtered.distance .- trt.distance[i]) 
                end
                # sort by difference to get the two closest
                sort!(filtered,[:diff])
                # skip if no control is available
                if nrow(filtered)==0
                    # add to unmatched
                    push!(unmatched,trt.id[i])
                    continue
                end
                # define quantities
                n_ctrl=min(nrow(filtered),ratio)
                # add to dataset
                pair=DataFrame(id=[trt.id[i];filtered.id[1:n_ctrl]],cluster=fill(i,n_ctrl+1))
                # push into dataset
                append!(matched,pair) 
            end
        end
    elseif distance == allowed_distances[3]
        # extract confounders from formula
        confounders = [t.sym for t in formula.rhs]
        # model matrix
        X=Matrix(df[!, confounders])
        # define mask
        mask=treatment_col .== true
        # treated and untreated
        X_ctrl = X[.!mask, :]
        # variance covariance matrix and regularization
        Σ = cov(X_ctrl, dims=1, corrected=false) + 1e-12*I
        # Σ positive definite
        if distance == "mahalanobis" && !isposdef(Σ)
            throw(ArgumentError("Covariance matrix not positive definite"))
        end
        # inverse
        Σ_inv = inv(Σ)
        # created filtered
        filtered=copy(df)
        if !replacement
            for r in 1:ratio
                for i in 1:nrow(trt)
                    # if nrow
                    # update from current filtered
                    current_mask = filtered[!, treatment] .== true
                    X = Matrix(filtered[!, confounders])
                    X_trt = X[current_mask, :]
                    X_untrt = X[.!current_mask, :]
                    # stop if untrt is empty
                    if size(X_untrt)[1]==0
                        @warn("Ran out of controls")
                        break
                    end
                    # calculate distance between current treated and untreated
                    mahalanobis_dists = calculate_mahalanobis(X_trt[i, :], X_untrt, Σ_inv)
                    idx_closest = argmin(mahalanobis_dists)
                    # extract IDs from CURRENT filtered
                    trt_patient_id = filtered[current_mask, :id][i]
                    ctrl_id = filtered[.!current_mask, :id][idx_closest]
                    # append
                    if r == 1
                        append!(matched, DataFrame(id=[trt_patient_id, ctrl_id], cluster=fill(i, 2)))
                    else
                        append!(matched, DataFrame(id=[ctrl_id], cluster=[i]))
                    end
                    # update filtered
                    filtered = filtered[filtered.id .≠ ctrl_id, :]
                end
            end
        else
            for i in 1:nrow(trt)
                # define matrixes
                X_trt = X[mask, :]
                X_untrt = X[.!mask, :]
                # calculate distance between current treated and untreated
                mahalanobis_dists = calculate_mahalanobis(X_trt[i, :], X_untrt, Σ_inv)
                closest_idxs = partialsortperm(mahalanobis_dists, 1:ratio)
                ctrl_ids = df[.!mask, :id][closest_idxs]
                # push to matched
                append!(matched, DataFrame(id=[trt.id[i]; ctrl_ids], cluster=fill(i, ratio+1)))
            end
        end
    end
    # warning if unmatched
    if length(unmatched)!=0
        @warn("Some treated are unmatched")
    end
    # if with replacement, delete cluster variable and keep unique rows
    if replacement
        # calculate weights
        ## number of controls per cluster
        cluster_size=combine(groupby(matched,:cluster),:cluster=>length=>:n_ctrl)
        ## create cluster weights (ALL units get 1/n_controls)
        matched.weight=1 ./ (cluster_size.n_ctrl[matched.cluster] .- 1)
        ## aggregate weights by ID (NO treated override)
        weight_table = combine(groupby(matched, :id), :weight => sum => :total_weight)
        ## merge to initial dataset
        leftjoin!(matched, weight_table, on=:id)
        # delete cluster column
        select!(matched,Not([:cluster,:weight]))
        # keep only unique rows
        matched_unique=unique(matched)
        # match with original dataset
        dfmatched=leftjoin(matched_unique,df, on=:id)
        # MatchIt normalization (AFTER aggregation)
        treated_mask = dfmatched[!, treatment] .== true
        norm_factor = mean(dfmatched[treated_mask, :total_weight])
        dfmatched[!, :total_weight] ./= norm_factor
        # sort by id
        sort!(dfmatched,:id)
    end
    # if replacement == false, delete those who do not have a complete set merge with original dataset
    if !replacement && discard==true
        dfmerged=leftjoin(matched,df, on=:id)
        counts=countmap(dfmerged.cluster)
        valid_clusters=[k for (k,v) in counts if v == ratio + 1]
        dfmatched=dfmerged[in.(dfmerged.cluster, Ref(valid_clusters)), :]
        sort!(dfmatched,:cluster)
    elseif !replacement && discard==false 
        dfmatched=leftjoin(matched,df, on=:id)
        sort!(dfmatched,:cluster)
    end    
    # sort based on cluster
    return PSMResult(dfmatched, unmatched, distance, ratio, replacement)
end
end