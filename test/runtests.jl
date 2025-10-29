using CausalTools, CSV, DataFrames, StatsModels
using Test

data_path = joinpath(@__DIR__, "..", "data", "df.csv")
df = CSV.read(data_path, DataFrame)

@testset "Inverse probability of treatment weighting tests" begin
    df_missing = deepcopy(df)
    allowmissing!(df_missing, :treatment)  # Make treatment column allow missing
    df_missing.treatment[1]=missing      # Now assign missing works
    iptw_obj = iptw(df, @formula(treatment~age+smoke+sex+cholesterol), truncate=[1,99], type="Stabilized")
    @test_throws AssertionError iptw(df_missing, @formula(treatment~age+smoke+sex+cholesterol), truncate=[1,99], type="Stabilized")
    @test typeof(iptw_obj) == CausalTools.IPTW.IPTWResult
    @test_throws AssertionError iptw(df, @formula(treatment ~ age + smoke + sex + cholesterol), truncate=[1, 99], type="Normalized")
    @test_throws AssertionError iptw(df, @formula(age ~ smoke + sex + cholesterol), truncate=[1, 99], type="Stabilized")
    @test_throws AssertionError iptw(df, @formula(treatment ~ smoke + sex + cholesterol), truncate=[180, 200], type="Stabilized")
end
