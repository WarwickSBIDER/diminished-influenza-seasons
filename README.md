# diminished-influenza-seasons
This repository houses code developed for the analysis presented in the scientific publication "Scenario modelling for diminished influenza seasons during 2020/2021 and 2021/2022 in England" by Edward M. Hill and Matt J. Keeling.

Model simulations are performed using the programming language Julia.
Julia makes use of environments, allowing bespoke package lists for separate projects. Documentation on working with environments and installing packages in the same state that is given by the project manifest: https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project

Please find below an explainer of the directory structure within this repository.  

## data
Comprised of the following subdirectories.

### AgeDistWithinRiskGrpWork
Age distribution data within at-risk and low risk groups.

### ContactData
Contact patterns by age group. Sourced from the study "Inferring the Structure of Social Contacts from Demographic Data in the Analysis of Infectious Diseases Spread" by L. Fumanelli et al. PLoS Comput Biol 8(9): e1002673. doi: 10.1371/journal.pcbi.1002673.

### DemographicData
Population size and mortality data (from ONS)

## HealthEconParamData
Summary data on the relative rate of occurrence (compared to GP consultations) for mortality, inpatient and outpatient attendances.

### RiskGrpPropnsData
By single year age groups, the proportion of the population of the given age that were in a seasonal influenza "at-risk' group. Season-by-season values spanning the 2000/01 - 2021/2022 influenza seasons.

### VaccEfficacy
Population-level vaccine efficacy estimates.

### VaccUptake
Daily vaccine uptake proportions for at-risk, low risk, entire population groups under differing vaccination scenarios.

## results
Directory containing simulation outputs and the plot script (household_bubble_impact_plots.m). The sub-directory 'matlab_packages' has files that are loaded by the plot script file

### 2020_2021_scenarios
Outputs associated with the scenario analysis conducted prior to the 2020/2021 influenza season.

**influenza_season_2020_2021_intervention_scenario_comparisons.jl**  
Script that produced the relative health episode occurrence MAT file, influenza_season_2020_2021_scenario_rel_vals.mat

**seasonal_influenza_2020_2021_scenarios_table_gen.m**  
Script used to compute the summary statistics shown in Tables XXX, YYY, ZZZ

**influenza_season_2020_2021_scenario_plots.m**  
Script used to produce Figures AAA, BBB, CCC.

### 2021_2022_scenarios
Outputs associated with the scenario analysis conducted prior to the 2021/2022 influenza season.

**influenza_season_2021_2022_intervention_scenario_comparisons.jl**  
Script that produced the relative health episode occurrence MAT file, influenza_season_2021_2022_scenario_rel_vals.mat

**seasonal_influenza_2021_2022_scenarios_table_gen.m**  
Script used to compute the summary statistics shown in Tables XXX, YYY, ZZZ

**seasonal_influenza_2021_2022_scenarios_table_gen.m**  
Script used to compute the summary statistics shown in Tables XXX, YYY, ZZZ

**influenza_season_2021_2022_scenario_plots.m**  
Script used to produce Figures AAA, BBB, CCC.

### matlab_packages
Files for packages that are used to produce and save the MATLAB generated figures.

## src

**RunSeasonalFluModelFMAlt_MultiSimn_2020_2021_scenarios.jl**
Script associated with running the scenario analysis conducted prior to the 2020/2021 influenza season. Generates outputs for the overall population (independent of risk status).

**RunSeasonalFluModelFMAlt_MultiSimn_2021_2022_scenarios.jl**
Script associated with running the scenario analysis conducted prior to the 2021/2022 influenza season. Generates outputs for the overall population (independent of risk status).

**RunSeasonalFluModelFMAlt_MultiSimn_RiskGrpVers_2020_2021_scenarios.jl**
Script associated with running the scenario analysis conducted prior to the 2020/2021 influenza season. Generates outputs for each risk status group.

**RunSeasonalFluModelFMAlt_MultiSimn_RiskGrpVers_2021_2022_scenarios.jl**
Script associated with running the scenario analysis conducted prior to the 2021/2022 influenza season. Generates outputs for each risk status group.

### model_supporting_functions
Functions to numerically solve ODEs and update exposure history array

### posterior_parameter_sets
Parameter sets used for parameterisation of the model are loaded from the files in this directory.
