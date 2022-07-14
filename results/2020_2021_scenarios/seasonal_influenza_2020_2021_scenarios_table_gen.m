%Purpose:
% Compute summary statistic values for seasonal influenza case severity
% outcomes showing relative impact of NPIs and expanded vaccination
% versus historical baseline
%--------------------------------------------------------------------------

clear variables

%% Load the relevant data file
load('influenza_season_2020_2021_scenario_rel_vals.mat','at_risk_rel_vals','low_risk_rel_vals','overall_rel_vals')

%% Set percentile values to be evaluated
percentile_vals = [50 2.5 97.5];

%% Generate summary statistic values for tables
all_popn_rel_vals_cell = compute_summ_stats_fn(overall_rel_vals,percentile_vals);
at_risk_rel_vals_cell = compute_summ_stats_fn(at_risk_rel_vals,percentile_vals);
low_risk_rel_vals_cell = compute_summ_stats_fn(low_risk_rel_vals,percentile_vals);

%% Compute relative value for each statistic and age group.
function output_cell = compute_summ_stats_fn(input_data,percentile_vals)

    %- SET GLOBAL VARIABLES & INTIALISE OUTPUT VARIABLES

    % Get number of statistics to be computed
    n_stats = numel(input_data);

    % Get number of age bins in use.
    % Given by row number of each three dimensional array
    n_age_bins = size(input_data{1},1);

    % Get number of non-baseline scenarios in use
    n_nonbaseline_scens = (size(input_data{1},3)) - 1;

    % Initialise output cell
    output_cell = cell(n_nonbaseline_scens,n_age_bins,n_stats);

    %- MAIN COMPUTATION LOOP

    % Iterate over each statistic.
    % Output computed values to the relevant row of storage cell.
    for stat_itr = 1:n_stats

        % Get data for this statistic
        stat_data = input_data{stat_itr};

        % Get number of scenarios in use. Number of slices of array
        n_scens = size(stat_data,3);

        % Iterate over each non-baseline scenario
        for scen_itr = 2:n_scens

            % Extract slice
            scen_data = stat_data(:,:,scen_itr);

            for age_grp_itr = 1:n_age_bins
                % Get age data for this statistic to be plotted
                config_data = scen_data(age_grp_itr,:);

                % Compute percentile values
                prctile_vals = prctile(config_data,percentile_vals);

                % Add to output cell
                output_cell{scen_itr-1,age_grp_itr,stat_itr} =...
                sprintf('%0.2f(%0.2f,%0.2f)',...
                            prctile_vals(1),... % Median
                            prctile_vals(2),... % Lower uncertainty bound
                            prctile_vals(3));   % Upper uncertainty bound
            end
        end
    end
end
