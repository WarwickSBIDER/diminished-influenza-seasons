#=
Purpose:
Take outputs from 2020/2021 influenza scenario runs and get
relative amounts of case severity counts
=#

#===========================
Set paths & load environment
===========================#

#Set paths
cd(dirname(@__FILE__))

#Load environment
using Pkg
Pkg.activate("../../")

#===========================
Load packages
===========================#
# Load libraries
using XLSX
using DelimitedFiles
using LinearAlgebra
using MAT


#===========================
Supporting functions for calculating proportion of single year age groups
experiencing given event severity
===========================#
# Function to construct linear ascertainment profile
function construct_ascertainment_profile(notch_ages,notch_vals)

    # Get linear piece lower and upper bounds
    AgeBandLowerBounds = notch_ages[1:end-1]
    AgeBandUpperBounds = notch_ages[2:end]

    # Inititalise storage vector for the profile
    n_ages = 101
    ascertainment_profile = zeros(n_ages)

    # Construct the profile
    for age_itr = 1:n_ages

        #Convert loop index to year of age
        YrOfAgeVal = age_itr - 1

        #Check upper bound, first <= than.
        AgeBandVectorIdx = findfirst(YrOfAgeVal .<= AgeBandUpperBounds)

        #Identify values at edge of bins, Entry idx and idx+1.
        LeftEdgeVal = notch_vals[AgeBandVectorIdx]
        RightEdgeVal = notch_vals[AgeBandVectorIdx+1]

        if AgeBandVectorIdx == length(AgeBandLowerBounds) #FOr last age category, get age years up to and including 100
            AgeBandSpan = 100 - AgeBandLowerBounds[AgeBandVectorIdx]
        else
          AgeBandSpan = AgeBandLowerBounds[AgeBandVectorIdx+1] - AgeBandLowerBounds[AgeBandVectorIdx]
        end

        IdxInAgeBand = (age_itr-1) - AgeBandLowerBounds[AgeBandVectorIdx] #Subtract 1 as ages begin from zero!
        ascertainment_profile[age_itr] = LeftEdgeVal + ((IdxInAgeBand/AgeBandSpan)*(RightEdgeVal - LeftEdgeVal))
    end

    return ascertainment_profile::Array{Float64,1}
end

# Fn to get proportion per year of age that are cases, hospitalsations and deaths under the simulated scenarios
# Output for single year of age. Can undergo further processing to be aggregated into requested age bands
function get_severity_counts(infected_propn_by_age,risk_status_val)

    # Set number of age groups in use
    n_ages = 101

    # Construct ascertainment profile
    notch_ages = [0,2,18,65,85,100]
    notch_vals = [0.2163,0.1422,0.2496,0.7022,0.5432,1]
    ascertainment_profile = construct_ascertainment_profile(notch_ages,notch_vals)

    # Get ascertained cases to benchmark other case severity against
    if risk_status_val == "low_risk"
        ascertained_by_age = 0.04.*ascertainment_profile
    elseif risk_status_val == "at_risk"
        AtRisk_scale = 1.5
        ascertained_by_age = 0.04.*ascertainment_profile.*AtRisk_scale
    else
        error("Invalid value for risk_status_val entered")
    end

    # Get 2020/2021 data
    cases_by_strain_2020_2021 = infected_propn_by_age[end,:,:]

    # Get ascertained cases per strain
    n_strains = 4
    ascertained_cases_by_age = zeros(n_ages,n_strains)
    for age_itr = 1:n_ages
        for strain_itr = 1:n_strains
            ascertained_cases_by_age[age_itr,strain_itr] = ascertained_by_age[age_itr].*cases_by_strain_2020_2021[age_itr,strain_itr]
        end
    end

    # Aggregate ascertained cases based on strain type
    ascertained_cases_type_A_2020_2021 = sum(ascertained_cases_by_age[:,1:2],dims=2)
    ascertained_cases_type_B_2020_2021 = sum(ascertained_cases_by_age[:,3:4],dims=2)
    ascertained_cases_2020_2021 = vec(sum(ascertained_cases_by_age,dims=2))

    # Inpatient admissions
    hospital_event_upper_age_bounds = [1,5,9,19,29,39,49,59,64,74,84,100]

    inpatient_low_risk_type_A = [0.7867,
0.131861,
0.056435,
0.058559,
0.090917,
0.074259,
0.058237,
0.064887,
0.060996,
0.112783,
0.306522,
1.264947,
]

    inpatient_at_risk_type_A = [4.160617,
0.884142,
0.341362,
0.183107,
0.185534,
0.220939,
0.271572,
0.334018,
0.363996,
0.743541,
1.005664,
1.469874
]

    inpatient_low_risk_type_B = [0.814172,
0.13299,
0.050888,
0,
0,
0,
0,
0,
0,
0,
0,
0,
]

    inpatient_at_risk_type_B = [3.462987,
0.853164,
0.323506,
0,
0,
0,
0,
0,
0,
0,
0,
0,
]

    # Populate the output vector
    inpatient_propn_of_age = zeros(n_ages)
    for age_itr = 1:n_ages
        yr_of_age = age_itr - 1

        bin_idx = findfirst(yr_of_age .<= hospital_event_upper_age_bounds)
        if risk_status_val == "low_risk"
            inpatient_propn_of_age[age_itr] =
                    (inpatient_low_risk_type_A[bin_idx].*ascertained_cases_type_A_2020_2021[age_itr]) +
                    (inpatient_low_risk_type_B[bin_idx].*ascertained_cases_type_B_2020_2021[age_itr])
        elseif risk_status_val == "at_risk"
            inpatient_propn_of_age[age_itr] =
                    (inpatient_at_risk_type_A[bin_idx].*ascertained_cases_type_A_2020_2021[age_itr]) +
                    (inpatient_at_risk_type_B[bin_idx].*ascertained_cases_type_B_2020_2021[age_itr])
        end
    end

    # In-hospital mortality
    #    No flu B mortality. Can sum over first two slices of array.

    # Set up relative event scalings
    inhospital_death_low_risk_type_A = [0.001982,
                                            0,
                                            0,
                                            0.000155,
                                            0.000112,
                                            8.81E-05,
                                            0.000458,
                                            0.000511,
                                            0.001452,
                                            0.002078,
                                            0.013218,
                                            0.130556]

    inhospital_death_at_risk_type_A = [0.202225,
                                        0.02782,
                                        0.009423,
                                        0.003846,
                                        0.003355,
                                        0.006459,
                                        0.015872,
                                        0.023081,
                                        0.032582,
                                        0.08246,
                                        0.136383,
                                        0.328687]

    # Populate the output vector
    in_hospital_death_propn_of_age = zeros(n_ages)
    for age_itr = 1:n_ages
        yr_of_age = age_itr - 1

        bin_idx = findfirst(yr_of_age .<= hospital_event_upper_age_bounds)
        if risk_status_val == "low_risk"
            in_hospital_death_propn_of_age[age_itr] =
                    inhospital_death_low_risk_type_A[bin_idx].*ascertained_cases_type_A_2020_2021[age_itr]
        elseif risk_status_val == "at_risk"
            in_hospital_death_propn_of_age[age_itr] =
                inhospital_death_at_risk_type_A[bin_idx].*ascertained_cases_type_A_2020_2021[age_itr]
        end
    end


    # Out of hospital mortality
        # Ages 50-64: 75% in hospital, 25% out of hospital
        # Ages 65-74: 65% in hospital, 35% out of hospital
        # Ages 75+: 50% in hospital, 50% out of hospital.
    out_of_hospital_death_propn_of_age = zeros(n_ages)
    for age_itr = 50:(n_ages-1)
        array_idx = age_itr + 1
        if age_itr < 65
            # Ages 50-64: 75% in hospital, 25% out of hospital
            out_of_hospital_death_propn_of_age[array_idx] = (25/75).*in_hospital_death_propn_of_age[array_idx]
        elseif age_itr < 75
            # Ages 65-74: 65% in hospital, 35% out of hospital
            out_of_hospital_death_propn_of_age[array_idx] = (35/65).*in_hospital_death_propn_of_age[array_idx]
        else
            # Ages 75+: 50% in hospital, 50% out of hospital.
            out_of_hospital_death_propn_of_age[array_idx] = 1. *in_hospital_death_propn_of_age[array_idx]
        end
    end

    return ascertained_cases_2020_2021::Array{Float64,1},
            inpatient_propn_of_age::Array{Float64,1},
            in_hospital_death_propn_of_age::Array{Float64,1},
            out_of_hospital_death_propn_of_age::Array{Float64,1}
end

# Produce the health event quantities for each risk group status
function get_event_propn_of_age_data(n_ages::Int64,n_replicates::Int64,
                                        infected_propn_by_age_atrisk,
                                        infected_propn_by_age_lowrisk)

    # Initialise the storage arrays
    at_risk_cases_propn_of_age = zeros(n_ages,n_replicates)
    at_risk_inpatient_propn_of_age = zeros(n_ages,n_replicates)
    at_risk_in_hospital_death_propn_of_age = zeros(n_ages,n_replicates)
    at_risk_out_of_hospital_death_propn_of_age = zeros(n_ages,n_replicates)

    low_risk_cases_propn_of_age = zeros(n_ages,n_replicates)
    low_risk_inpatient_propn_of_age = zeros(n_ages,n_replicates)
    low_risk_in_hospital_death_propn_of_age = zeros(n_ages,n_replicates)
    low_risk_out_of_hospital_death_propn_of_age = zeros(n_ages,n_replicates)

    # Populate arrays with data from each replicate
    for replicate_itr = 1:n_replicates
        # At risk
        risk_status_val = "at_risk"
        at_risk_cases_propn_of_age[:,replicate_itr],
        at_risk_inpatient_propn_of_age[:,replicate_itr],
        at_risk_in_hospital_death_propn_of_age[:,replicate_itr],
        at_risk_out_of_hospital_death_propn_of_age[:,replicate_itr] = get_severity_counts(infected_propn_by_age_atrisk[replicate_itr],
                                                                                    risk_status_val)
        # Low risk
        risk_status_val = "low_risk"
        low_risk_cases_propn_of_age[:,replicate_itr],
        low_risk_inpatient_propn_of_age[:,replicate_itr],
        low_risk_in_hospital_death_propn_of_age[:,replicate_itr],
        low_risk_out_of_hospital_death_propn_of_age[:,replicate_itr] = get_severity_counts(infected_propn_by_age_lowrisk[replicate_itr],
                                                                    risk_status_val)
    end

    # Get estimates for risk groups combined
    overall_cases_propn_of_age = overall_popn_propn_of_age_estimates(low_risk_cases_propn_of_age,
                                                                                at_risk_cases_propn_of_age)
    overall_inpatient_propn_of_age = overall_popn_propn_of_age_estimates(low_risk_inpatient_propn_of_age,
                                                                            at_risk_inpatient_propn_of_age)
    overall_in_hospital_death_propn_of_age = overall_popn_propn_of_age_estimates(low_risk_in_hospital_death_propn_of_age,
                                                                                at_risk_in_hospital_death_propn_of_age)
    overall_out_of_hospital_death_propn_of_age = overall_popn_propn_of_age_estimates(low_risk_out_of_hospital_death_propn_of_age,
                                                                            at_risk_out_of_hospital_death_propn_of_age)


    return at_risk_cases_propn_of_age,
            at_risk_inpatient_propn_of_age,
            at_risk_in_hospital_death_propn_of_age,
            at_risk_out_of_hospital_death_propn_of_age,
            low_risk_cases_propn_of_age,
            low_risk_inpatient_propn_of_age,
            low_risk_in_hospital_death_propn_of_age,
            low_risk_out_of_hospital_death_propn_of_age,
            overall_cases_propn_of_age,
            overall_inpatient_propn_of_age,
            overall_in_hospital_death_propn_of_age,
            overall_out_of_hospital_death_propn_of_age
end

#===========================
Supporting functions for recombining risk group estimates into overall population
===========================#
# Take inputs from low risk and at risk and output estimates for total population
# at that single year of age
function overall_popn_propn_of_age_estimates(low_risk_data::Array{Float64,2},
                                                at_risk_data::Array{Float64,2})

    # Load proportion of each year of age that is at-risk
    atrisk_propn_per_age = XLSX.readdata("../../data/RiskGrpPropnsData/AtRiskPropnSpreadsheets/ByYrOfAge_RiskGrpPropnPerSeason_ModifiedAges2to4_EMH.xlsx", "Sheet1", "B21:CX21")
    atrisk_propn_per_age = convert(Array{Float64,1},atrisk_propn_per_age[:])

    # Combine fractions from atrisk and lowrisk, with appropriate scaling
    n_ages = size(low_risk_data,1)
    n_replicates = size(low_risk_data,2)
    overall_propn_per_age = zeros(n_ages,n_replicates)
    for replicate_itr = 1:n_replicates
        for age_itr = 1:n_ages
            overall_propn_per_age[age_itr,replicate_itr] = (at_risk_data[age_itr,replicate_itr]*atrisk_propn_per_age[age_itr]) +
                                                                (low_risk_data[age_itr,replicate_itr]*(1 - atrisk_propn_per_age[age_itr]))
        end
    end
    return overall_propn_per_age::Array{Float64,2}
end

#===========================
Supporting function for performing age aggregation
===========================#
# Take in the data. Aggregate accroding to age_agg_bin_upper_bounds values
# As final age bin, calculate estimates for all ages
function aggregate_data_into_age_classes(age_agg_bin_upper_bounds::Array{Int64,1},
                                        popn_by_yr_of_age::Array{Float64,1},
                                        input_data_tuple::Array{Array{Array{Float64,2},1},1},
                                        n_replicates::Int64
                                        )

    # Get number of age bins to be used
    n_age_bins = length(age_agg_bin_upper_bounds) + 1 # Add one for all ages estimates

    # Get number of statistics to be iterated over
    n_stats = length(input_data_tuple)

    # Get index to be accessed, based on year of ages that are start and end point of each age bin
    start_index = [1;age_agg_bin_upper_bounds[1:end-1].+2]
    end_index = age_agg_bin_upper_bounds.+1

    # Initialise the output tuple
    output_data_tuple = Array{Array{Array{Float64,2},1},1}(undef,n_stats)

    # Populate the output tuple for each statistic
    for stat_itr = 1:n_stats

        # Extract the relevant input data for this iteration
        current_itr_input_data = input_data_tuple[stat_itr]

        # Get number of scenarios to be compared
        n_scens = length(current_itr_input_data)

        # Initialise the array for current statistic
        output_data_tuple[stat_itr] = Array{Array{Float64,2},1}(undef,n_scens)

        for scen_itr = 1:n_scens

            # Assign data for current scenario to variable
            current_scen_data = current_itr_input_data[scen_itr]

            # Initialise output array for this scenario & statistic
            output_data_tuple[stat_itr][scen_itr] = zeros(n_age_bins,n_replicates)

            # Compute statistic for the age bin
            for age_bin_itr = 1:n_age_bins

                if age_bin_itr == n_age_bins
                    this_age_bin_start_idx = 1
                    this_age_bin_end_idx = end_index[end]
                else
                    # Get index values of population array to be used for this age bin
                    this_age_bin_start_idx = start_index[age_bin_itr]
                    this_age_bin_end_idx = end_index[age_bin_itr]
                end

                # Get relative size of each single year of age as propn of the current age bin
                single_yr_propn_of_age_bin = popn_by_yr_of_age[this_age_bin_start_idx:this_age_bin_end_idx]./sum(popn_by_yr_of_age[this_age_bin_start_idx:this_age_bin_end_idx])

                # Add up values attributed to each single year of age, scaled by single_yr_propn_of_age_bin
                for replicate_itr = 1:n_replicates
                    output_data_tuple[stat_itr][scen_itr][age_bin_itr,replicate_itr] = sum(current_scen_data[this_age_bin_start_idx:this_age_bin_end_idx,replicate_itr].*single_yr_propn_of_age_bin)
                end
            end
        end
    end

    # Return the vector of vector of arrays
    return output_data_tuple::Array{Array{Array{Float64,2},1},1}
end

#===========================
Supporting functions for computing relative counts (compared to baseline strategy)
===========================#
# Compare outputs from the different scenarios and get values relative to baseline
function compare_scenarios_to_baseline(age_agg_data_tuple::Array{Array{Array{Float64,2},1},1})

    # Get number of statistics
    n_stats = length(age_agg_data_tuple)

    # Initialise the output variable. Relative values per statistic
    # Within each statistic, a three dimensional array
    # Row per age bin, column per replicate, slice per scenario compared to baseline
    rel_vals = Array{Array{Float64,3},1}(undef,n_stats)

    # For each statistic, get number of scenarios.
    # First scenario is baseline. Store that. Iterate over others and compare to baseline
    for stat_itr = 1:n_stats
        current_stat_data = age_agg_data_tuple[stat_itr]
        n_scens = length(current_stat_data)

        baseline_data = current_stat_data[1]
        n_age_bins::Int64 = size(baseline_data,1)
        n_replicates::Int64 = size(baseline_data,2)
        current_stat_rel_val = zeros(n_age_bins,n_replicates,n_scens)
        for scen_itr = 1:n_scens
            current_stat_rel_val[:,:,scen_itr] = current_stat_data[scen_itr]./baseline_data
        end

        # Assign array to output variable
        rel_vals[stat_itr] = copy(current_stat_rel_val)
    end

    return rel_vals::Array{Array{Float64,3},1}
end

#===========================
Main script begins here
===========================#

#===========================
Specify variables consistent across scenarios
===========================#
n_ages = 101
n_replicates = 100

#===========================
Load the baseline data
===========================#
# Low risk
MATfile_baseline_lowrisk = matopen("FMAlt_SimnOutputFiles_Julia/FMAlt_SimnJuliaV1ModelRun_#Covid19_baseline_LowRisk.mat")
SimnDataAllPopn_baseline_lowrisk = read(MATfile_baseline_lowrisk,"SimnData")
close(MATfile_baseline_lowrisk)

infected_propn_by_age_baseline_lowrisk = SimnDataAllPopn_baseline_lowrisk[:,1]
popn_dist_baseline_lowrisk = SimnDataAllPopn_baseline_lowrisk[:,2]

# At risk
MATfile_baseline_atrisk = matopen("FMAlt_SimnOutputFiles_Julia/FMAlt_SimnJuliaV1ModelRun_#Covid19_baseline_AtRisk.mat")
SimnDataAllPopn_baseline_atrisk = read(MATfile_baseline_atrisk,"SimnData")
close(MATfile_baseline_atrisk)

infected_propn_by_age_baseline_atrisk = SimnDataAllPopn_baseline_atrisk[:,1]
popn_dist_baseline_atrisk = SimnDataAllPopn_baseline_atrisk[:,2]

#===========================
Load the expanded vaccine programme only scenario data
===========================#
# Low risk
MATfile_expandvacc_lowrisk = matopen("FMAlt_SimnOutputFiles_Julia/FMAlt_SimnJuliaV1ModelRun_#Covid19_vacc_only_LowRisk.mat")
SimnDataAllPopn_expandvacc_lowrisk = read(MATfile_expandvacc_lowrisk,"SimnData")
close(MATfile_expandvacc_lowrisk)

infected_propn_by_age_expandvacc_lowrisk = SimnDataAllPopn_expandvacc_lowrisk[:,1]
popn_dist_expandvacc_lowrisk = SimnDataAllPopn_expandvacc_lowrisk[:,2]

# At risk
MATfile_expandvacc_atrisk = matopen("FMAlt_SimnOutputFiles_Julia/FMAlt_SimnJuliaV1ModelRun_#Covid19_vacc_only_AtRisk.mat")
SimnDataAllPopn_expandvacc_atrisk = read(MATfile_expandvacc_atrisk,"SimnData")
close(MATfile_expandvacc_atrisk)

infected_propn_by_age_expandvacc_atrisk = SimnDataAllPopn_expandvacc_atrisk[:,1]
popn_dist_expandvacc_atrisk = SimnDataAllPopn_expandvacc_atrisk[:,2]

#===========================
Load the NPIs only scenario data (with vaccination matching what has occurred historically)
===========================#
# Low risk
MATfile_NPIs_lowrisk = matopen("FMAlt_SimnOutputFiles_Julia/FMAlt_SimnJuliaV1ModelRun_#Covid19_NPIs_only_LowRisk.mat")
SimnDataAllPopn_NPIs_lowrisk = read(MATfile_NPIs_lowrisk,"SimnData")
close(MATfile_NPIs_lowrisk)

infected_propn_by_age_NPIs_lowrisk = SimnDataAllPopn_NPIs_lowrisk[:,1]
popn_dist_NPIs_lowrisk = SimnDataAllPopn_NPIs_lowrisk[:,2]

# At risk
MATfile_NPIs_atrisk = matopen("FMAlt_SimnOutputFiles_Julia/FMAlt_SimnJuliaV1ModelRun_#Covid19_NPIs_only_AtRisk.mat")
SimnDataAllPopn_NPIs_atrisk = read(MATfile_NPIs_atrisk,"SimnData")
close(MATfile_NPIs_atrisk)

infected_propn_by_age_NPIs_atrisk = SimnDataAllPopn_NPIs_atrisk[:,1]
popn_dist_NPIs_atrisk = SimnDataAllPopn_NPIs_atrisk[:,2]

#===========================
Load the combined NPIs and vaccine expansion scenario data
===========================#
# Low risk
MATfile_expandvacc_and_NPIs_lowrisk = matopen("FMAlt_SimnOutputFiles_Julia/FMAlt_SimnJuliaV1ModelRun_#Covid19_NPIs_and_vacc_LowRisk.mat")
SimnDataAllPopn_expandvacc_and_NPIs_lowrisk = read(MATfile_expandvacc_and_NPIs_lowrisk,"SimnData")
close(MATfile_expandvacc_and_NPIs_lowrisk)

infected_propn_by_age_expandvacc_and_NPIs_lowrisk = SimnDataAllPopn_expandvacc_and_NPIs_lowrisk[:,1]
popn_dist_expandvacc_and_NPIs_lowrisk = SimnDataAllPopn_expandvacc_and_NPIs_lowrisk[:,2]

# At risk
MATfile_expandvacc_and_NPIs_atrisk = matopen("FMAlt_SimnOutputFiles_Julia/FMAlt_SimnJuliaV1ModelRun_#Covid19_NPIs_and_vacc_AtRisk.mat")
SimnDataAllPopn_expandvacc_and_NPIs_atrisk = read(MATfile_expandvacc_and_NPIs_atrisk,"SimnData")
close(MATfile_expandvacc_and_NPIs_atrisk)

infected_propn_by_age_expandvacc_and_NPIs_atrisk = SimnDataAllPopn_expandvacc_and_NPIs_atrisk[:,1]
popn_dist_expandvacc_and_NPIs_atrisk = SimnDataAllPopn_expandvacc_and_NPIs_atrisk[:,2]

#===========================
Process the data
===========================#
# Baseline
at_risk_baseline_cases_propn_of_age,
        at_risk_baseline_inpatient_propn_of_age,
        at_risk_baseline_in_hospital_death_propn_of_age,
        at_risk_baseline_out_of_hospital_death_propn_of_age,
        low_risk_baseline_cases_propn_of_age,
        low_risk_baseline_inpatient_propn_of_age,
        low_risk_baseline_in_hospital_death_propn_of_age,
        low_risk_baseline_out_of_hospital_death_propn_of_age,
        overall_baseline_cases_propn_of_age,
        overall_baseline_inpatient_propn_of_age,
        overall_baseline_in_hospital_death_propn_of_age,
        overall_baseline_out_of_hospital_death_propn_of_age = get_event_propn_of_age_data(n_ages,n_replicates,
                                                                                    infected_propn_by_age_baseline_atrisk,
                                                                                    infected_propn_by_age_baseline_lowrisk)

# Expanded vaccination only
at_risk_expandvacc_cases_propn_of_age,
        at_risk_expandvacc_inpatient_propn_of_age,
        at_risk_expandvacc_in_hospital_death_propn_of_age,
        at_risk_expandvacc_out_of_hospital_death_propn_of_age,
        low_risk_expandvacc_cases_propn_of_age,
        low_risk_expandvacc_inpatient_propn_of_age,
        low_risk_expandvacc_in_hospital_death_propn_of_age,
        low_risk_expandvacc_out_of_hospital_death_propn_of_age,
        overall_expandvacc_cases_propn_of_age,
        overall_expandvacc_inpatient_propn_of_age,
        overall_expandvacc_in_hospital_death_propn_of_age,
        overall_expandvacc_out_of_hospital_death_propn_of_age = get_event_propn_of_age_data(n_ages,n_replicates,
                                                                                    infected_propn_by_age_expandvacc_atrisk,
                                                                                    infected_propn_by_age_expandvacc_lowrisk)


# NPIs only
at_risk_NPIs_cases_propn_of_age,
        at_risk_NPIs_inpatient_propn_of_age,
        at_risk_NPIs_in_hospital_death_propn_of_age,
        at_risk_NPIs_out_of_hospital_death_propn_of_age,
        low_risk_NPIs_cases_propn_of_age,
        low_risk_NPIs_inpatient_propn_of_age,
        low_risk_NPIs_in_hospital_death_propn_of_age,
        low_risk_NPIs_out_of_hospital_death_propn_of_age,
        overall_NPIs_cases_propn_of_age,
        overall_NPIs_inpatient_propn_of_age,
        overall_NPIs_in_hospital_death_propn_of_age,
        overall_NPIs_out_of_hospital_death_propn_of_age = get_event_propn_of_age_data(n_ages,n_replicates,
                                                                                    infected_propn_by_age_NPIs_atrisk,
                                                                                    infected_propn_by_age_NPIs_lowrisk)

# Expanded vaccination & NPIs
at_risk_expandvacc_and_NPIs_cases_propn_of_age,
        at_risk_expandvacc_and_NPIs_inpatient_propn_of_age,
        at_risk_expandvacc_and_NPIs_in_hospital_death_propn_of_age,
        at_risk_expandvacc_and_NPIs_out_of_hospital_death_propn_of_age,
        low_risk_expandvacc_and_NPIs_cases_propn_of_age,
        low_risk_expandvacc_and_NPIs_inpatient_propn_of_age,
        low_risk_expandvacc_and_NPIs_in_hospital_death_propn_of_age,
        low_risk_expandvacc_and_NPIs_out_of_hospital_death_propn_of_age,
        overall_expandvacc_and_NPIs_cases_propn_of_age,
        overall_expandvacc_and_NPIs_inpatient_propn_of_age,
        overall_expandvacc_and_NPIs_in_hospital_death_propn_of_age,
        overall_expandvacc_and_NPIs_out_of_hospital_death_propn_of_age = get_event_propn_of_age_data(n_ages,n_replicates,
                                                                                    infected_propn_by_age_expandvacc_and_NPIs_atrisk,
                                                                                    infected_propn_by_age_expandvacc_and_NPIs_lowrisk)

#===========================
Put outputs to be compared into vector of arrays
===========================#
# Put data to be compared into vector of arrays
low_risk_cases_propn_of_age = [low_risk_baseline_cases_propn_of_age,low_risk_expandvacc_cases_propn_of_age,
                                    low_risk_NPIs_cases_propn_of_age, low_risk_expandvacc_and_NPIs_cases_propn_of_age]
low_risk_inpatient_propn_of_age = [low_risk_baseline_inpatient_propn_of_age,low_risk_expandvacc_inpatient_propn_of_age,
                                        low_risk_NPIs_inpatient_propn_of_age, low_risk_expandvacc_and_NPIs_inpatient_propn_of_age]
low_risk_in_hospital_death_propn_of_age = [low_risk_baseline_in_hospital_death_propn_of_age,low_risk_expandvacc_in_hospital_death_propn_of_age,
                                                low_risk_NPIs_in_hospital_death_propn_of_age, low_risk_expandvacc_and_NPIs_in_hospital_death_propn_of_age]
low_risk_out_of_hospital_death_propn_of_age = [low_risk_baseline_out_of_hospital_death_propn_of_age,low_risk_expandvacc_out_of_hospital_death_propn_of_age,
                                                low_risk_NPIs_out_of_hospital_death_propn_of_age, low_risk_expandvacc_and_NPIs_out_of_hospital_death_propn_of_age]


low_risk_data_tuple = [low_risk_cases_propn_of_age,low_risk_inpatient_propn_of_age,
                            low_risk_in_hospital_death_propn_of_age,low_risk_out_of_hospital_death_propn_of_age]

at_risk_cases_propn_of_age = [at_risk_baseline_cases_propn_of_age,at_risk_expandvacc_cases_propn_of_age,
                                    at_risk_NPIs_cases_propn_of_age, at_risk_expandvacc_and_NPIs_cases_propn_of_age]
at_risk_inpatient_propn_of_age = [at_risk_baseline_inpatient_propn_of_age,at_risk_expandvacc_inpatient_propn_of_age,
                                        at_risk_NPIs_inpatient_propn_of_age, at_risk_expandvacc_and_NPIs_inpatient_propn_of_age]
at_risk_in_hospital_death_propn_of_age = [at_risk_baseline_in_hospital_death_propn_of_age,at_risk_expandvacc_in_hospital_death_propn_of_age,
                                                at_risk_NPIs_in_hospital_death_propn_of_age, at_risk_expandvacc_and_NPIs_in_hospital_death_propn_of_age]
at_risk_out_of_hospital_death_propn_of_age= [at_risk_baseline_out_of_hospital_death_propn_of_age,at_risk_expandvacc_out_of_hospital_death_propn_of_age,
                                                at_risk_NPIs_out_of_hospital_death_propn_of_age, at_risk_expandvacc_and_NPIs_out_of_hospital_death_propn_of_age]


at_risk_data_tuple = [at_risk_cases_propn_of_age,at_risk_inpatient_propn_of_age,
                            at_risk_in_hospital_death_propn_of_age,at_risk_out_of_hospital_death_propn_of_age]

overall_cases_propn_of_age = [overall_baseline_cases_propn_of_age,overall_expandvacc_cases_propn_of_age,
                                    overall_NPIs_cases_propn_of_age, overall_expandvacc_and_NPIs_cases_propn_of_age]
overall_inpatient_propn_of_age = [overall_baseline_inpatient_propn_of_age,overall_expandvacc_inpatient_propn_of_age,
                                        overall_NPIs_inpatient_propn_of_age, overall_expandvacc_and_NPIs_inpatient_propn_of_age]
overall_in_hospital_death_propn_of_age = [overall_baseline_in_hospital_death_propn_of_age,overall_expandvacc_in_hospital_death_propn_of_age,
                                            overall_NPIs_in_hospital_death_propn_of_age, overall_expandvacc_and_NPIs_in_hospital_death_propn_of_age]
overall_out_of_hospital_death_propn_of_age= [overall_baseline_out_of_hospital_death_propn_of_age,overall_expandvacc_out_of_hospital_death_propn_of_age,
                                            overall_NPIs_out_of_hospital_death_propn_of_age, overall_expandvacc_and_NPIs_out_of_hospital_death_propn_of_age]

overall_data_tuple = [overall_cases_propn_of_age,overall_inpatient_propn_of_age,
                            overall_in_hospital_death_propn_of_age,overall_out_of_hospital_death_propn_of_age]
#===========================
Put outputs to be compared into vector of arrays
Put single year of age data through age aggregation
Compare outputs between scenarios
===========================#

# Specify the upper limit of each age band to be used
age_agg_bin_upper_bounds = [1,17,49,64,79,100]

# Load the single year of age population data
low_risk_popn_by_yr_of_age = readdlm("../../data/DemographicData/ONSPopnDistEngland20102019_0to100_LowRisk.txt",',')[end,:]
at_risk_popn_by_yr_of_age = readdlm("../../data/DemographicData/ONSPopnDistEngland20102019_0to100_AtRisk.txt",',')[end,:]
overall_popn_by_yr_of_age = readdlm("../../data/DemographicData/ONSPopnDistEngland20102019_0to100.txt",',')[end,:]

# Low risk
age_agg_low_risk_data_tuple = aggregate_data_into_age_classes(age_agg_bin_upper_bounds,
                                                            low_risk_popn_by_yr_of_age,
                                                            low_risk_data_tuple,
                                                            n_replicates)
# At risk
age_agg_at_risk_data_tuple =  aggregate_data_into_age_classes(age_agg_bin_upper_bounds,
                                    at_risk_popn_by_yr_of_age,
                                    at_risk_data_tuple,
                                    n_replicates)
# Overall
age_agg_overall_data_tuple =  aggregate_data_into_age_classes(age_agg_bin_upper_bounds,
                                    overall_popn_by_yr_of_age,
                                    overall_data_tuple,
                                    n_replicates)

#===========================
Compare outputs between scenarios
===========================#
# Across the scenarios, take values relative to the baseline scenario
low_risk_rel_vals = compare_scenarios_to_baseline(age_agg_low_risk_data_tuple)
at_risk_rel_vals = compare_scenarios_to_baseline(age_agg_at_risk_data_tuple)
overall_rel_vals = compare_scenarios_to_baseline(age_agg_overall_data_tuple)

#===========================
Save outputs to file
===========================#
# Have MAT with each variable saved to the same file
file = matopen("influenza_season_2020_2021_scenario_rel_vals_copy.mat", "w")
write(file, "low_risk_rel_vals", low_risk_rel_vals)
write(file, "at_risk_rel_vals", at_risk_rel_vals)
write(file, "overall_rel_vals", overall_rel_vals)
close(file)
