using CausalTools, Test, Random, DataFrames, CSV
using StatsModels: @formula

data_path = joinpath(@__DIR__, "..", "data", "df.csv")
df = CSV.read(data_path, DataFrame)

@testset "Inverse probability of treatment weighting tests" begin
    df_missing = deepcopy(df)
    allowmissing!(df_missing, :treatment) 
    df_missing.treatment[1]=missing    
    iptw_obj = iptw(df, @formula(treatment~age+smoke+sex+cholesterol), truncate=[1,99], type="Stabilized")
    @test_throws AssertionError iptw(df_missing, @formula(treatment~age+smoke+sex+cholesterol), truncate=[1,99], type="Stabilized")
    @test typeof(iptw_obj) == CausalTools.IPTW.IPTWResult
    @test_throws AssertionError iptw(df, @formula(treatment ~ age + smoke + sex + cholesterol), truncate=[1, 99], type="Normalized")
    @test_throws AssertionError iptw(df, @formula(age ~ smoke + sex + cholesterol), truncate=[1, 99], type="Stabilized")
    @test_throws AssertionError iptw(df, @formula(treatment ~ smoke + sex + cholesterol), truncate=[180, 200], type="Stabilized")
end

@testset "Propensity score matching tests" begin
    # test invalid inputs (mirroring iptw structure)
    df_missing = deepcopy(df)
    allowmissing!(df_missing, :treatment) 
    df_missing.treatment[1] = missing
    @test_throws ArgumentError psm(df_missing, @formula(treatment ~ age + smoke + sex + cholesterol))
    @test_throws ArgumentError psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), ratio=0)
    @test_throws ArgumentError psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), ratio=1.5)
    @test_throws ArgumentError psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), replacement="true")
    @test_throws ArgumentError psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), distance="invalid")
    @test_throws ArgumentError psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), caliper=-0.1)
    @test_throws ArgumentError psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), order="invalid")
    @test_throws ArgumentError psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), discard="true")
    @test_throws ArgumentError psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), seed=1.5)
    # test invalid formula (wrong lhs)
    @test_throws ArgumentError psm(df, @formula(age ~ smoke + sex + cholesterol))
    # test missing covariates (mahalanobis)
    @test_throws ArgumentError psm(df, @formula(treatment ~ age + smoke + sex + cholesterol + missing_var), distance="mahalanobis")
    # test no treated/no controls
    df_no_trt = deepcopy(df)
    df_no_trt.treatment .= false
    @test_throws ArgumentError psm(df_no_trt, @formula(treatment ~ age + smoke + sex + cholesterol))
    df_no_ctrl = deepcopy(df)
    df_no_ctrl.treatment .= true
    @test_throws ArgumentError psm(df_no_ctrl, @formula(treatment ~ age + smoke + sex + cholesterol))
    # test non-binary treatment
    df_nonbinary = deepcopy(df)
    @test_throws ArgumentError psm(df_nonbinary, @formula(blood_pressure ~ age + smoke + sex + cholesterol))
    # test valid psm object creation (default)
    psm_obj = psm(df, @formula(treatment ~ age + smoke + sex + cholesterol))
    @test typeof(psm_obj) == CausalTools.PSM.PSMResult
    @test psm_obj.distance == "probability"
    @test psm_obj.ratio == 1
    @test psm_obj.replacement == false
    @test nrow(psm_obj.dataset) > 0
    @test eltype(psm_obj.unmatched) == Int
    @test hasproperty(psm_obj.dataset, :id)
    @test hasproperty(psm_obj.dataset, :cluster)
    # test different configurations
    psm_logit = psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), distance="logit")
    @test psm_logit.distance == "logit"
    #psm_mahal = psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), distance="mahalanobis")
    #@test psm_mahal.distance == "mahalanobis"
    psm_ratio2 = psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), ratio=2)
    @test psm_ratio2.ratio == 2
    psm_replace = psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), replacement=true)
    @test psm_replace.replacement == true
    @test hasproperty(psm_replace.dataset, :total_weight) 
    @test all(psm_replace.dataset.total_weight .> 0)
    psm_discard = psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), ratio=2, discard=true)
    @test psm_discard.ratio == 2
    # test ordering effects (different results)
    psm_largest = psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), order="largest")
    psm_smallest = psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), order="smallest")
    @test psm_largest.dataset != psm_smallest.dataset
    # test random seed reproducibility
    Random.seed!(123)
    psm_rand1 = psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), order="random", seed=123)
    Random.seed!(123)
    psm_rand2 = psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), order="random", seed=123)
    @test psm_rand1.dataset == psm_rand2.dataset
    # test caliper reduces matches
    psm_nocaliper = psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), caliper=Inf)
    psm_caliper = psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), caliper=0.1)
    @test nrow(psm_nocaliper.dataset) >= nrow(psm_caliper.dataset)
    # test discard removes incomplete clusters
    psm_keep_partial = psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), ratio=5, discard=false)
    psm_discard_complete = psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), ratio=5, discard=true)
    @test nrow(psm_keep_partial.dataset) >= nrow(psm_discard_complete.dataset)
end
