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

function psm(data::DataFrame, formula::FormulaTerm; distance::Nothing="Logistic",method::Nothing="Greedy", caliper::Union{Nothing,Integer,Float64}=nothing,ratio::Integer=1,replacement::Bool=false,link::Link=LogitLink())
    # copy data
    df=copy(data)
    # define id
    df[!,:id]=1:nrow(df)
    # extract treatment from left side function
    treatment=formula.lhs.sym
    # make sure caliper is between 0 and 1

    # make sure treatment is binary
    treatment_col=df[!,treatment]
    @assert !any(ismissing, treatment_col) "Treatment column contains missing values."
    @assert all(v -> v in (0,1), treatment_col) "Treatment must be coded as 0/1."
    # define distance to match on, whether mahalabonis or probability
    if distance=="Logistics"
        trt_model=glm(formula,df,Binomial(),link)
        df.distance=predict(trt_model)
        df.expit=log.((df.probability) ./ (1 .- df.probability))
    elseif  distance=="Mahalabonis"
    end
    # divide initial dataset into treated and untreated
    untrt=data[data[!,exposure].==false,:]
    trt=data[data[!,exposure] .==true,:]
    # define index for untrt and trt to delete from dataset those that went throught the loop
    untrt[!, :index] = 1:nrow(untrt)
    trt[!, :index] = 1:nrow(trt)
    # sort treated patients
    sort!(trt, :distance,rev=true)
    # initialize matched dataset
    clustered=DataFrame(cluster=Int[],id=Int[])
    # initialize unmatched vector 
    unmatched=Int[]; lowerthanratio=Int[] 
    # rows to loop through
    rows=size(trt)[1]
    # perform matching
     for j in 1:ratio
        for i in 1:rows
            # break if untrt dataset is empty
            if nrow(untrt) == 0
                println("Warning: You ran out of controls, some treated patients will not be matched")
                append!(unmatched,trt[i:end, :id])
                break  
            end
            index_trt=trt[i,:index]
            untrti=copy(untrt)
            untrti.diff=abs.(trt.distance[i] .- untrti.distance)
            # if caliper is not equal to nothing, then filter
            if caliper==nothing
                
            end 
            # 
            if nrow(untrti)!=0 && nrow(untrti)<ratio
                println("Warning: Not enough matches for treated unit $i")
                push!(lowerthanratio,i)
                append!(clustered, DataFrame(cluster=[i],id=[trt.id[i]]))
                continue
            end
            if nrow(untrti) ==0
                println("Warning: No matches for treated unit $i")
                push!(unmatched,i)
                append!(clustered, DataFrame(cluster=[i],id=[trt.id[i]]))
                continue
            end
            # arrange based on difference
            sort!(untrti,[:diff])
            indexes=untrti[1,:index]
            # create new df with the new observations and the cluster identificator
            if j ==1
                append!(clustered, DataFrame(cluster=[i],id=[trt.id[i]]))
            end
            append!(clustered, DataFrame(cluster=[i],id=[untrti.id[1]]))
            if replace==false
                untrt=filter(row -> !(row.index in indexes), untrt)  
            end 
        end
    end
    # merge clustered dataset to initial DataFrame
    matched=innerjoin(clustered,data,on=:id,makeunique=true)
    return matched, unmatched, lowerthanratio
end