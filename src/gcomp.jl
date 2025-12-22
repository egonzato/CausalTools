using GLM, Statistics, StatsBase, LSurvival, LinearAlgebra

# predict survival probabilities up to time T
function predictsurv(model,dataset)
    # extract covariates
    covars=[t.sym for t in model.formula.rhs.terms]
    # extract bh from model
    bh=model.bh[:,1]
    # extract coefficients
    ## coefficients
    coeffs=coef(model)
    # create dataset for predictions
    X=dataset[:,covars]
    # linear predictors
    lp=[dot(dataset[i, covars], coeffs) for i in 1:nrow(dataset)] 
    #lp=lp .- mean(lp)
    # probabilities
    S_t=exp.(-bh .* exp.(lp'))
    # return matrix
    return S_t, model.bh[:,4], lp
end
# g-computation


function gcomp(dataset,formula;treatment=:treatment,bootstrap=false,family="gaussian",nbootstrap=1000,contrast="difference")
    # copy dataset
    df=copy(dataset)
    # create treated and untreated
    df_trt=copy(df); df_trt[!, treatment] .= true 
    df_ctrl=copy(df); df_ctrl[!, treatment] .= false
    # based on family argument, define model 
    if family=="gaussian"
        # model
        model=lm(formula,df)
        # predictions
        pred_trt=predict(model, df_trt)
        pred_ctrl=predict(model, df_ctrl)
        # calculate ate
        ate=mean(pred_trt-pred_ctrl)
        # bootstrap if needed
        bootate=[]
        se=NaN
        if bootstrap
            bootate=[]
            for i in 1:nbootstrap
                bootdf=df[sample(1:nrow(df),nrow(df),replace=true),:]
                # model
                boot_model=lm(formula,bootdf)
                # create treated and untreated
                df_trt=copy(bootdf); df_trt[!, treatment] .= true 
                df_ctrl=copy(bootdf); df_ctrl[!, treatment] .= false
                # predict 
                pred_trt=predict(boot_model, df_trt)
                pred_ctrl=predict(boot_model, df_ctrl)
                # calculate ate
                push!(bootate,mean(pred_trt-pred_ctrl))
            end
            # calculate std of ate distribution
            se=std(bootate)
        end
    elseif family=="binomial"
        # model
        model=glm(formula,df,Binomial())
        # predictions
        pred_trt=predict(model, df_trt)
        pred_ctrl=predict(model, df_ctrl)
        # calculate ate
        if contrast=="difference"
            ate=mean(pred_trt-pred_ctrl)
        else 
            ate=mean(pred_trt)/mean(pred_ctrl)
        end
        # bootstrap if needed
        bootate=[]
        se=NaN
        if bootstrap
            bootate=[]
            for i in 1:nbootstrap
                bootdf=df[sample(1:nrow(df),nrow(df),replace=true),:]
                # model
                boot_model=glm(formula,bootdf,Binomial())
                # create treated and untreated
                df_trt=copy(bootdf); df_trt[!, treatment] .= true 
                df_ctrl=copy(bootdf); df_ctrl[!, treatment] .= false
                # predict 
                pred_trt=predict(boot_model, df_trt)
                pred_ctrl=predict(boot_model, df_ctrl)
                # calculate ate
                if contrast=="difference"
                    current_bootate=mean(pred_trt-pred_ctrl)
                else 
                    current_bootate=mean(pred_trt)/mean(pred_ctrl)
                end
                # calculate ate
                push!(bootate,current_bootate)
            end
            # calculate std of ate distribution
            se=std(bootate)
        end
    else
        # model
        model=coxph(formula,df,ties = "efron")
        # predict on treated
        pred, time, lps=predictsurv(model,df_trt)
        # keep last time for treated
        pred_trt=pred[end,:]
        # do the same but for untreated
        # predict on treated
        pred, time, lps=predictsurv(model,df_ctrl)
        # keep last time for treated
        pred_ctrl=pred[end,:]
        # ate
        if contrast=="difference"
            ate=mean(pred_trt .- pred_ctrl)
        else 
            #ate=mean(pred_trt)/mean(pred_ctrl)
            ate = log.((.- log.(mean(pred_trt))) ./ (.- log.(mean(pred_ctrl))))
        end
        # bootstrap if needed
        bootate=[]
        se=NaN
        if bootstrap
            bootate=[]
            for i in 1:nbootstrap
                bootdf=df[sample(1:nrow(df),nrow(df),replace=true),:]
                # model
                boot_model=coxph(formula,bootdf,ties = "efron")
                # create treated and untreated
                df_trt=copy(bootdf); df_trt[!, treatment] .= true 
                df_ctrl=copy(bootdf); df_ctrl[!, treatment] .= false
                # predict 
                # predict on treated
                pred, time, lps=predictsurv(boot_model,df_trt)
                # keep last time for treated
                pred_trt=pred[end,:]
                # do the same but for untreated
                # predict on treated
                pred, time, lps=predictsurv(boot_model,df_ctrl)
                # keep last time for treated
                pred_ctrl=pred[end,:]
                # calculate ate
                if contrast=="difference"
                    current_bootate=mean(pred_trt .- pred_ctrl)
                else 
                    #ate=mean(pred_trt)/mean(pred_ctrl)
                    current_bootate = log.((.- log.(mean(pred_trt))) ./ (.- log.(mean(pred_ctrl))))
                end
                # calculate ate
                push!(bootate,current_bootate)
            end
            # calculate std of ate distribution
            se=std(bootate)
        end
    end
    # what should the function give back?
    return ate, model, bootstrap, se, bootate 
end