%Purpose:
% Produce violin plots showing relative impact of a lack of influenza season
% on the epidemiological outcomes in the subsequent influenza season

% Counterfactual: Transmission not impacted in prior influenza season.

% Date: June 2021
% Modified: July 2022
%--------------------------------------------------------------------------

clear variables


%% ADD PATH DEPENDENCIES
addpath('../matlab_packages/Violinplot-Matlab')
addpath('../matlab_packages/export_fig')

%% Load the relevant data file
load('influenza_season_2021_2022_scenarios_rel_vals.mat','at_risk_rel_vals','low_risk_rel_vals','overall_rel_vals')

%% Set up plotting variables

% Set x-positions data will be plotted at
% Not plotting baseline scenario, so first row of x_plot_pos corresponds to scenario 2 etc
x_plot_pos = [1:7;10:16;19:25];

% Set colour values for each age bracket
colour_vals = [0, 0.4470, 0.7410;
                0.8500, 0.3250, 0.0980;
                0.9290, 0.6940, 0.1250;
                0.4940, 0.1840, 0.5560;
                0.4660, 0.6740, 0.1880;
                0.3010, 0.7450, 0.9330;
                0.5, 0.5, 0.5];

% Set up legend labels
legend_label = {'0-1yrs', '2-17yrs','18-49yrs','50-64yrs','65-84yrs','85+yrs','All ages'};

% Set x-axis labels
xaxis_label = 'With suppressed influenza season';
xticks_vals = 4;
xticks_labels = {''};

%xaxis_label = 'Scenarios';
%xticks_vals = [4 13 22];
%xticks_labels = {'Expanded vacc.','NPIs','NPIs & expanded vacc.'};

%Set axes limits
x_limits = [0 8];
y_limits = [0.0 2];

% Set titles
title_vec = {'Symptomatic (consult GP)',...
                'Inpatient admissions',...
                'In-hospital mortality',...
                'Out of hospital mortality'};

% Set plot fontsize
plot_fontsize = 20;

%% Call the multi-panel figure function for each risk group

% Low risk
low_risk_fig_title = 'Low-risk';
low_risk_save_filename = 'low_risk_comparison';
plot_violin_multi_panels(low_risk_rel_vals,...
                            x_plot_pos,...
                            colour_vals,...
                            legend_label,...
                            xaxis_label,...
                            xticks_vals,...
                            xticks_labels,...
                            x_limits,...
                            y_limits,...
                            title_vec,...
                            plot_fontsize,...
                            low_risk_fig_title,...
                            low_risk_save_filename)
% At risk
at_risk_fig_title = 'At-risk';
at_risk_save_filename = 'at_risk_comparison';
plot_violin_multi_panels(at_risk_rel_vals,...
                            x_plot_pos,...
                            colour_vals,...
                            legend_label,...
                            xaxis_label,...
                            xticks_vals,...
                            xticks_labels,...
                            x_limits,...
                            y_limits,...
                            title_vec,...
                            plot_fontsize,...
                            at_risk_fig_title,...
                            at_risk_save_filename)

% Overall
overall_fig_title = 'Independent of risk status';
overall_save_filename = 'overall_comparison';
plot_violin_multi_panels(overall_rel_vals,...
                            x_plot_pos,...
                            colour_vals,...
                            legend_label,...
                            xaxis_label,...
                            xticks_vals,...
                            xticks_labels,...
                            x_limits,...
                            y_limits,...
                            title_vec,...
                            plot_fontsize,...
                            overall_fig_title,...
                            overall_save_filename)

%% Call the single panel figure function

% Low risk
low_risk_save_filename_vec = {'single_panel_figures/low_risk_comparison_symptomatic';...
                             'single_panel_figures/low_risk_comparison_inpatient_admissions';...
                              'single_panel_figures/low_risk_comparison_inhospital_mortality';...
                              'single_panel_figures/low_risk_comparison_out_of_hospital_mortality'};
plot_violin_single_panel(low_risk_rel_vals,...
                            x_plot_pos,...
                            colour_vals,...
                            legend_label,...
                            xaxis_label,...
                            xticks_vals,...
                            xticks_labels,...
                            x_limits,...
                            y_limits,...
                            title_vec,...
                            plot_fontsize,...
                            low_risk_save_filename_vec)

% At risk
at_risk_save_filename_vec = {'single_panel_figures/at_risk_comparison_symptomatic';...
                             'single_panel_figures/at_risk_comparison_inpatient_admissions';...
                              'single_panel_figures/at_risk_comparison_inhospital_mortality';...
                              'single_panel_figures/at_risk_comparison_out_of_hospital_mortality'};
plot_violin_single_panel(at_risk_rel_vals,...
                            x_plot_pos,...
                            colour_vals,...
                            legend_label,...
                            xaxis_label,...
                            xticks_vals,...
                            xticks_labels,...
                            x_limits,...
                            y_limits,...
                            title_vec,...
                            plot_fontsize,...
                            at_risk_save_filename_vec)

% Overall
overall_save_filename_vec = {'single_panel_figures/overall_comparison_symptomatic';...
                             'single_panel_figures/overall_comparison_inpatient_admissions';...
                              'single_panel_figures/overall_comparison_inhospital_mortality';...
                              'single_panel_figures/overall_comparison_out_of_hospital_mortality'};
plot_violin_single_panel(overall_rel_vals,...
                            x_plot_pos,...
                            colour_vals,...
                            legend_label,...
                            xaxis_label,...
                            xticks_vals,...
                            xticks_labels,...
                            x_limits,...
                            y_limits,...
                            title_vec,...
                            plot_fontsize,...
                            overall_save_filename_vec)

%% Plotting function for multi-panel figures

function plot_violin_multi_panels(input_data,...
                            x_plot_pos,...
                            colour_vals,...
                            legend_label,...
                            xaxis_label,...
                            xticks_vals,...
                            xticks_labels,...
                            x_limits,...
                            y_limits,...
                            title_vec,...
                            plot_fontsize,...
                            figure_title,...
                            save_filename)

  % Initialise figures
    position = [10, 10, 2.5*550, 2.5*450];
    set(0, 'DefaultFigurePosition', position);
    %set(0, 'DefaultFigurePosition', position,'Units','points');
    fig = figure();
    clf;
    set(fig,'Color', [1 1 1])
    hold on

  % Get number of statistics to be plotted
  % Separate panel for each
  n_stats = numel(input_data);

  % Get number of age bins in use. Given by row number of each three
  % dimensional array
  n_age_bins = size(input_data{1},1);

  % Produce 2x2 set of panels. One panel per statistic
  for stat_itr = 1:n_stats
    subplot(2,2,stat_itr)
    hold on

    % Get data for this statistic
    stat_data = input_data{stat_itr};

    % Get number of scenarios in use. Number of slices of array
    n_scens = size(stat_data,3);

    % Pick out relevant data for this statistic
    stat_data = input_data{stat_itr};
    for scen_itr = 2:n_scens
        % Iterate over each non-baseline scenario. Plot data for each age
        % bin
        scen_data = stat_data(:,:,scen_itr);
        for age_bin_itr = 1:n_age_bins
            % Get age data for this statistic to be plotted
            config_data = scen_data(age_bin_itr,:);

            % Get colour for the plot
            colour_vec = colour_vals(age_bin_itr,:);

            %Create violin plot for each age band
            if (stat_itr ~= 4) || (age_bin_itr > 3) % Conditions to avoid plotting NaN for young ages that have zero out of hospital mortality
                % Otherwise, create violin plot
                 violins = Violin(config_data, x_plot_pos(scen_itr-1,age_bin_itr));
                                             % Subtract one as not plotting baseline.

                 % Set violin plot properties
                 violins.ViolinColor = colour_vec;
                 violins.ViolinAlpha = 1.0; %shading transparency

                 %Set violin plot region properties
                 violins.EdgeColor = [1 1 1];

                 %Set median marker properties
                 violins.MedianColor = [1 1 1];
                 violins.MedianPlot.Marker = 's';
                 violins.MedianPlot.MarkerEdgeColor = [0 0 0];
                 violins.MedianPlot.LineWidth = 1;

                 %Set whisker line properties
                 violins.BoxColor = [0 0 0];
            end
        end
    end

    % Add line at y=1
    plot([0 x_limits(end)],[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',1.5);

    % Set axes limits
    xlim(x_limits)
    ylim(y_limits)

    % Add x-axis label, if required
    if stat_itr > 2
        xlabel(xaxis_label)
    end

    % Add xtick labels
    xticks(xticks_vals)
    xticklabels(xticks_labels)

    % Add y-axis labels
    if (stat_itr == 1) || (stat_itr == 3)
        ylabel('Relative to no COVID-19 counterfactual')
    end

    % Add title
    title(title_vec{stat_itr})

    % Add legend
    if stat_itr == 2
        % Initialise patch vector and populate
        H = gobjects(n_age_bins,1);
        for age_bin_itr = 1:n_age_bins
          H(age_bin_itr) = fill(NaN,NaN,colour_vals(age_bin_itr,:),'DisplayName',legend_label{age_bin_itr});
        end

        % Construct the legend
        legend(H,...
                      'LineWidth',1.5,...
                      'FontSize',plot_fontsize,...
                      'Orientation','horizontal',...
                      'Position',[0.33890909090909 0.520444444444445 0.567636363636363 0.0248888888888888])
    end

    %Specify general axis properties
      set(gca,'FontSize',plot_fontsize)
      set(gca,'LineWidth',1)
      box on
  end

  % Add overall figure title
  sgtitle(figure_title,'FontWeight','bold','FontSize',plot_fontsize+4)

  % Save figure to file
  export_fig(save_filename,'-pdf','-r1200')
end


%% Plotting function for single panel figures

function plot_violin_single_panel(input_data,...
                                    x_plot_pos,...
                                    colour_vals,...
                                    legend_label,...
                                    xaxis_label,...
                                    xticks_vals,...
                                    xticks_labels,...
                                    x_limits,...
                                    y_limits,...
                                    title_vec,...
                                    plot_fontsize,...
                                    save_filename_vec)


  % Get number of statistics to be plotted
  % Separate panel for each
  n_stats = numel(input_data);

  % Get number of age bins in use. Given by row number of each three
  % dimensional array
  n_age_bins = size(input_data{1},1);

  % Produce a figure per statistic
  for stat_itr = 1:n_stats

    % Initialise figure
    position = [10, 10, 1.2*550, 1.2*450];
    set(0, 'DefaultFigurePosition', position);
    %set(0, 'DefaultFigurePosition', position,'Units','points');
    fig = figure();
    clf;
    set(fig,'Color', [1 1 1])
    hold on

    % Get data for this statistic
    stat_data = input_data{stat_itr};

    % Get number of scenarios in use. Number of slices of array
    n_scens = size(stat_data,3);

    % Pick out relevant data for this statistic
    stat_data = input_data{stat_itr};
    for scen_itr = 2:n_scens
        % Iterate over each non-baseline scenario. Plot data for each age
        % bin
        scen_data = stat_data(:,:,scen_itr);
        for age_bin_itr = 1:n_age_bins
            % Get age data for this statistic to be plotted
            config_data = scen_data(age_bin_itr,:);

            % Get colour for the plot
            colour_vec = colour_vals(age_bin_itr,:);

            %Create violin plot for each age band
            if (stat_itr ~= 4) || (age_bin_itr > 3) % Conditions to avoid plotting NaN for young ages that have zero out of hospital mortality
                % Otherwise, create violin plot
                 violins = Violin(config_data, x_plot_pos(scen_itr-1,age_bin_itr));
                                             % Subtract one as not plotting baseline.

                 % Set violin plot properties
                 violins.ViolinColor = colour_vec;
                 violins.ViolinAlpha = 1.0; %shading transparency

                 %Set violin plot region properties
                 violins.EdgeColor = [1 1 1];

                 %Set median marker properties
                 violins.MedianColor = [1 1 1];
                 violins.MedianPlot.Marker = 's';
                 violins.MedianPlot.MarkerEdgeColor = [0 0 0];
                 violins.MedianPlot.LineWidth = 1;

                 %Set whisker line properties
                 violins.BoxColor = [0 0 0];
            end
        end
    end

    % Add line at y=1
    plot([0 x_limits(end)],[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',1.5);

    % Set axes limits
    xlim(x_limits)
    ylim(y_limits)

    % Add x-axis label, if required
    xlabel(xaxis_label)

    % Add xtick labels
    xticks(xticks_vals)
    xticklabels(xticks_labels)

    % Add y-axis labels
    ylabel('Relative to no COVID-19 counterfactual')

    % Add title
    title(title_vec{stat_itr})

    % Add legend
    if stat_itr == 2
        % Initialise patch vector and populate
        H = gobjects(n_age_bins,1);
        for age_bin_itr = 1:n_age_bins
            H(age_bin_itr) = fill(NaN,NaN,colour_vals(age_bin_itr,:),'DisplayName',legend_label{age_bin_itr});
        end

        % Construct the legend
        legend(H,...
            'LineWidth',1.5,...
            'FontSize',plot_fontsize,...
            'Orientation','vertical',...
            'Location','southeast');
        %'Position',[0.33890909090909 0.520444444444445 0.567636363636363 0.0248888888888888])
    end

    %Specify general axis properties
    set(gca,'FontSize',plot_fontsize)
    set(gca,'LineWidth',1)
    box on

    % Save figure to file
    export_fig(save_filename_vec{stat_itr},'-pdf','-r1200')
  end
end
