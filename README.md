## Roadmap / To do

### 0. Package setup

- [X] Rename `data/df.CSV` → `data/df.csv` (Linux/CI is case-sensitive)
- [X] Clean `Project.toml`: remove `Revise`, `Documenter` (and unused `Colors`, `IterTools`, `Distributions`, `CSV`) from `[deps]`; move test-only packages to `[extras]`
- [X] Fix `[compat]`: standard libraries pinned to 1.11 conflict with `julia = "1.6.7"` (either raise minimum Julia to 1.10 or relax the entries)
- [X] Remove `docs/build/` from git and add it to `.gitignore`
- [ ] Get CI green
- [ ] Consider moving plotting into a package extension, so `Plots` isn't a hard dependency

### 1. Matching (`psm`)

**Bug fixes**

- [ ] Symbol vs String comparison with `names(df)`: Mahalanobis always throws, and `id` is added to the user's dataset
- [ ] Never modify the input DataFrame
- [ ] Check inputs (binary treatment, missing values) *before* fitting the GLM
- [ ] Ratio > 1 without replacement: controls can be added to clusters whose treated unit was never matched; `unmatched` gets duplicates
- [ ] "Ran out of controls" should stop all rounds, not just the inner loop
- [ ] Replacement branch: use the `treatment` variable instead of the hard-coded `:treatment`; agree on weight column names (`weight` / `nweight` vs `total_weight` in tests)
- [X] Re-enable the Mahalanobis test

**Features**

- [ ] Caliper for Mahalanobis distance (raw caliper and/or within a PS caliper)
- [ ] `link` argument for the propensity model (logit / probit / cloglog), reusing `check_link`
- [ ] Support categorical variables and interactions in Mahalanobis (use `modelmatrix`)
- [ ] Ratio ≥ 1 (fixed) and, optionally, variable ratio
- [ ] Estimand argument (ATT / ATC)
- [ ] **Optimal matching**: Hungarian algorithm for 1:1, min-cost flow for 1:k
- [ ] Exact matching on selected variables (optional)
- [ ] Performance: avoid copying/sorting the DataFrame for each treated unit

### 2. Inverse probability of treatment weighting (`iptw`)

- [ ] Estimands: ATE (done), ATT, ATC, overlap weights
- [ ] Truncation: decide between truncating weights and trimming the PS; fix the docstring ("untreated units")
- [ ] Effective sample size
- [ ] Effect estimation: weighted outcome model with robust (sandwich) and bootstrap SE
- [ ] Optional: categorical or continuous treatment

### 3. Diagnostics (shared by matching and IPTW)

- [ ] Standardized mean differences before and after matching/weighting
- [ ] Variance ratios
- [ ] Love plot
- [ ] Propensity score overlap plot
- [ ] `show` methods for `PSMResult` / `IPTWResult` summarizing what was done

### 4. G-computation (`gcomp`) — optional

- [ ] Include and export it from `CausalTools.jl`
- [ ] Return a `GCompResult` struct instead of a tuple
- [ ] Factor out the repeated bootstrap code
- [ ] Percentile and normal bootstrap confidence intervals
- [ ] ATE / ATT
- [ ] Survival outcomes: add `LSurvival` as a dependency (or postpone), and check the Cox contrast

### 5. Tests and documentation

- [ ] Tests that check results, not only errors (e.g. compare against R's `MatchIt` / `WeightIt` on `df.csv`)
- [ ] Docstrings for all exported functions
- [ ] Documenter site with one worked example per method
- [ ] Update the README with `psm` and `gcomp` examples
