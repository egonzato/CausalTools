using CausalTools, CSV, DataFrames, StatsModels
using Test

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
    # Test invalid inputs (mirroring IPTW structure)
    df_missing = deepcopy(df)
    allowmissing!(df_missing, :treatment)
    df_missing.treatment[1] = missing
    
    @test_throws AssertionError psm(df_missing, @formula(treatment ~ age + smoke + sex + cholesterol))
    @test_throws AssertionError psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), ratio=0.5)
    @test_throws AssertionError psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), ratio=2, replacement="true")
    @test_throws AssertionError psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), ratio=2, distance="invalid")
    @test_throws AssertionError psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), ratio=2, caliper=-0.1)
    @test_throws AssertionError psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), ratio=2, order="invalid")
    @test_throws AssertionError psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), ratio=2, discard="true")
    
    # Test invalid formula (non-binary treatment)
    df_nonbinary = deepcopy(df)
    df_nonbinary.treatment[1] = 2
    @test_throws ArgumentError psm(df_nonbinary, @formula(treatment ~ age + smoke + sex + cholesterol))
    
    # Test valid PSM object creation
    psm_obj = psm(df, @formula(treatment ~ age + smoke + sex + cholesterol))
    @test typeof(psm_obj) == CausalTools.PSM.PSMResult
    @test psm_obj.distance == "probability"
    @test psm_obj.ratio == 1
    @test psm_obj.replacement == false
    @test nrow(psm_obj.dataset) > 0
    @test eltype(psm_obj.unmatched) == Int
    
    # Test different configurations
    psm_mahal = psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), distance="mahalanobis")
    @test psm_mahal.distance == "mahalanobis"
    
    psm_ratio2 = psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), ratio=2)
    @test psm_ratio2.ratio == 2
    
    psm_replace = psm(df, @formula(treatment ~ age + smoke + sex + cholesterol), replacement=true)
    @test psm_replace.replacement == true
    

    @test hasproperty(psm_obj.dataset, :id)
    @test hasproperty(psm_obj.dataset, :cluster) || psm_obj.replacement  
end
