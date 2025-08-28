%----------------------------------------------------------------------
% Post-processing and visualization of KGW spectral analysis results
%----------------------------------------------------------------------
clear all, close all, clc

% Figure settings (common to all plots)
set(groot,'defaultAxesFontSize',20) % figures font size
set(groot,'DefaultTextFontSize',20) % figures font size
set(groot, 'defaultAxesTickLabelInterpreter','latex'); % latex labels
set(groot, 'defaultLegendInterpreter','latex'); % latex labels
set(groot, 'defaultTextInterpreter','latex'); % latex labels

% Configuration parameters
config.n = 3;                  % Excitation index (must match solver)
config.nEigval = 8;             % Number of eigenvalues in convergence plot
config.nEigmod = 4;             % Number of eigenmodes in v(r), w(r) plots (4 for n=0,1,2)
config.whichEigmod = 'first';    % Plot 'last' eigenmodes or 'first' eigenmodess
config.cmpEigmod = 'on';        % Compare eigenmodes with background stationary states ('on' or 'off')
config.saveFigures = true;      % Flag to control figure saving


% Input/Output parameters
io.inputFolder_res = "Results";
io.outputFolder_fig = "Figures";
io.inputFile_res = "KGWspectral" + sprintf('_n%03d', config.n) + ".mat";


% Figure output files
io.outputFile_fig_specIdx = "KGWspectral" + sprintf('_n%03d', config.n) + "_spectrum_Idx.png";
io.outputFile_fig_specConv = "KGWspectral" + sprintf('_n%03d', config.n) + "_spectrum_Conv.png";
io.outputFile_fig_specRelDiff = "KGWspectral" + sprintf('_n%03d', config.n) + "_relativeDiff.png";
io.outputFile_fig_eigv = "KGWspectral" + sprintf('_n%03d', config.n) + "_eigenmodes_v_"+config.whichEigmod+"_"+config.cmpEigmod+".png";
io.outputFile_fig_eigw = "KGWspectral" + sprintf('_n%03d', config.n) + "_eigenmodes_w_"+config.whichEigmod+"_"+config.cmpEigmod+".png";

% Graphics parameters
config.colorVec_default = [
    [0      0.4470 0.7410];
    [0.8500 0.3250 0.0980];
    [0.9290 0.6940 0.1250];
    [0.4940 0.1840 0.5560];
    [0.4660 0.6740 0.1880];
    [0.3010 0.7450 0.9330];
    [0.6350 0.0780 0.1840]
];


%--- ANALYSIS AT FIXED BACKGROUND n -------------------------------------

% Load Results at fixed n
[params, results_Nvary] = loadResults(io);

if config.n>3
    % 1. Plot spectrum vs index
    plotSpectrumIndex(results_Nvary, io, config, config.saveFigures);
end

% 2.a. Plot spectrum convergence vs N
plotSpectrumConvergence(results_Nvary, io, config, config.saveFigures);

% 2.b. Plot relative differences for convergence analysis
plotConvergenceAnalysis(results_Nvary, config, io, config.saveFigures);

% 3. Plot eigenmodes V(r)
plotEigenmodes(results_Nvary, config, io, config.colorVec_default, 'V', config.saveFigures);

% 4. Plot eigenmodes W(r)
plotEigenmodes(results_Nvary, config, io, config.colorVec_default, 'W', config.saveFigures);

% 5. Print Stability Analysis
printStabilityAnalysis(results_Nvary,config)
%------------------------------------------------------------------------





%%

%----------------------------------------------------------------------
%--- PLOTTING FUNCTIONS ----------------------------------------------
%----------------------------------------------------------------------

function [params, results_Nvary] = loadResults(io)
    % Load KGW spectral results 

    % Create output folder
    createSubfolder("./" + io.outputFolder_fig);

    % Load results
    results_file = io.inputFolder_res + "/" + io.inputFile_res;
    if ~isfile(results_file)
        error('Results file not found: %s\n', results_file);
    end
    
    data = load(results_file);
    results_Nvary = data.results_Nvary;
    params = data.params;
    
    fprintf('Loaded results for n=%d with N values: [%s]\n', ...
        params.n, sprintf('%d ', params.Nvec));

end

function plotSpectrumIndex(results_Nvary, io, config, saveFig)
    % Plot eigenvalue spectrum  as a function of the index (unstable only)
    % Use N max as a reference

    figure()
    set(gcf,'Position',[100 100 560*1 420])
    hold on

    % Reference eigenvalues
    ref_omega2 = results_Nvary(end).omega2 ;
    % Retrieve stable indices, sort descending
    [desc_omega2,idx_desc] = sort(ref_omega2,'descend');
    unstable_idx = results_Nvary(end).stability_analysis.unstable_idx(idx_desc);
    unstable_idx_progressive = cumsum(unstable_idx);
    % Define xdata, ydata
    xdata = unstable_idx_progressive(unstable_idx) -1;  %-1 since k starts from 0
    ydata = desc_omega2(unstable_idx);


    % Fit tests
    % 1) exponential: -y = b exp(ax)  
    p = polyfit(xdata, log(-ydata), 1);
    plot(xdata, -exp(xdata).^p(1) .*exp(p(2)), ...
        "Color",[config.colorVec_default(1,:), 0.5],'LineWidth',1.5,'DisplayName','exponential')
    % 2) parabolic: y = ax^2 +bx +c
    p = polyfit(xdata, ydata, 2);
    plot(xdata, p(1)*xdata.^2 + p(2)*xdata+p(3),...
        "Color",[config.colorVec_default(2,:), 0.5],'LineWidth',1.5,'DisplayName','parabolic')
    % 3) cubic: y = ax^3 +bx^2 +cx +d
    p = polyfit(xdata, ydata, 3);
    plot(xdata, p(1)*xdata.^3+ p(2)*xdata.^2 + p(3)*xdata+p(4),...
        "Color",[config.colorVec_default(3,:), 0.5],'LineWidth',1.5,'DisplayName','cubic')
    % % 4) power law: -y = b x^a
    % p = polyfit(log(xdata(2:0)), log(-ydata(2:end)), 1);
    % plot(xdata(2:end), -xdata(2:end).^p(1) .*exp(p(2)),...
    %     "Color",[config.colorVec_default(4,:), 0.5],'LineWidth',1.5,'DisplayName','power')

    % Add data plot
    plot(xdata, ydata, '.',"MarkerFaceColor",config.colorVec_default(1,:), 'MarkerSize', 8, 'DisplayName', 'Unstable');

    % Format Plot
    xlabel('$k$');
    ylabel('$\omega^2$');
    leg = legend('Location', 'southwest');
    set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
    % Axis Limits
    % [~,ref_idx] = min(abs( ref_omega2 - 4*abs(min(ref_omega2)) ));
    % xlim([0, stable_idx_progressive(ref_idx)]);
    % ylim([2*min(ref_omega2), 4*min(abs(ref_omega2))]);

    % Save Figure
    if saveFig
        saveas(gcf, io.outputFolder_fig + "/" + io.outputFile_fig_specIdx);
        fprintf('  Saved: %s\n', io.outputFile_fig_specIdx);
    end
end

function plotSpectrumConvergence(results_Nvary, io, config, saveFig)
    % Plot eigenvalue spectrum as function of discretization N (unstable only)
    
    figure();
    set(gcf,'Position',[100 100 560*1 420]) 
    hold on;

    % Extract results across different N to plot all at once  
    stuck_unstable_idx = arrayfun(@(x) x.stability_analysis.unstable_idx , results_Nvary, 'UniformOutput', false); 
    stuck_omega2 = {results_Nvary.omega2}; 
    stuck_N = arrayfun(@(x) x.N * ones(size(x.omega2)) , results_Nvary, 'UniformOutput', false);
    % Concatenate into single vector 
    stuck_unstable_idx = vertcat(stuck_unstable_idx{:}); 
    stuck_omega2 = vertcat(stuck_omega2{:}); 
    stuck_N = vertcat(stuck_N{:});

    % Plot all unstable modes at once
    if sum(stuck_unstable_idx)>0
        plot(stuck_omega2(stuck_unstable_idx), stuck_N(stuck_unstable_idx), '.',"MarkerFaceColor",config.colorVec_default(1,:), 'MarkerSize', 8, 'DisplayName', 'Unstable');
    end
    
    % Format plot
    xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlabel('$\omega^2$');
    ylabel('$N$');
    
    % Set axis limits
    stuck_omega2 = results_Nvary(1).omega2;
    xlim([2*min(stuck_omega2), 4*abs(min(stuck_omega2))]);
    ylim([2*results_Nvary(1).N - results_Nvary(2).N, 2*results_Nvary(end).N - results_Nvary(end-1).N]);
    
    grid on;
    leg = legend('Location', 'northeast');
    set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
    
    if saveFig
        saveas(gcf, io.outputFolder_fig + "/" + io.outputFile_fig_specConv);
        fprintf('  Saved: %s\n', io.outputFile_fig_specConv);
    end
end

function plotConvergenceAnalysis(results_Nvary, config, io, saveFig)
    % Plot relative differences between consecutive N values (unstable only)
    
    % Retrieve N vector
    Nvec = [results_Nvary.N]';
    % Calculate relative differences
    N0 = Nvec(1);
    diff_matrix = zeros(N0, length(Nvec) - 1);    
    for i = 2:length(results_Nvary)
        omega2_prev = results_Nvary(i-1).omega2(1:N0);
        omega2_curr = results_Nvary(i).omega2(1:N0);
        omega2_ref = results_Nvary(1).omega2(1:N0);
        
        diff_matrix(:, i-1) = (omega2_curr - omega2_prev) ./ omega2_ref;
    end
    
    % Plot
    figure();
    set(gcf,'Position',[100 100 560*1 420]) 
    hold on;
    
    % nr of plots
    nr_unstable_modes = results_Nvary(1).stability_analysis.nr_unstable_modes;
    nr_plots = min(config.nEigval, nr_unstable_modes);
    
    % Plot differences
    for i = 1:nr_plots
        % Blue with increasing transparency for increasing indices 
        alpha = 1.0 - (i-1) / nr_plots;
        color = config.colorVec_default(1,:);
        
        % Plot differences
        scatter(Nvec(2:end), diff_matrix(i+nr_unstable_modes-nr_plots, :), '+', ... %or, for first modes: diff_matrix(i, :)
            'MarkerEdgeColor', color, 'MarkerEdgeAlpha', alpha, ...
             'DisplayName', sprintf('$\\omega^2_{%d}$', i+nr_unstable_modes-nr_plots -1), 'LineWidth', 1.5);   %-1 since k starts from 0    %or, for first modes: i
    end
    % Plot parabolic fit
    fprintf("\n=== Convergence analysis, n="+string(config.n)+" ===\n")
    for i = 1:nr_plots
        % Blue with increasing transparency for increasing indices 
        alpha = 1.0 - (i-1) / nr_plots;
        color = config.colorVec_default(1,:);

        % Parabolic fit: y = ax^2 +bx +c
        p = polyfit(Nvec(2:end), diff_matrix(i+nr_unstable_modes-nr_plots, :), 2);  %or, for first modes: diff_matrix(i, :)
        plot(Nvec(2:end), p(1)*Nvec(2:end).^2 + p(2)*Nvec(2:end)+p(3), '-', 'Color', [color(:); alpha], 'HandleVisibility','off');%'DisplayName', sprintf("$\\Delta = %+.2eN^2%+.2eN%+.2e$",p))
   
        % Print to command window
        fprintf(sprintf("mode "+string(i-1)+",\tDelta = %+.2e N^2  %+.2e N  %+.2e\n",p))  %-1 since k starts from 0
    end
    
    xlabel('$N$');
    ylabel('$\Delta $');%(N_j) = \frac{\omega^2(N_j)-\omega^2(N_{j-1})}{\omega^2(N_0)}$');
    leg = legend('Location', 'southeast','NumColumns',2);
    set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
    grid on;
    
    if saveFig
        saveas(gcf, io.outputFolder_fig + "/" + io.outputFile_fig_specRelDiff);
        fprintf('\n  Saved: %s\n', io.outputFile_fig_specRelDiff);
    end
end

function plotEigenmodes(results_Nvary, config, io, colorVec_default, component, saveFig)
    % Plot eigenmode components (V or W) as function of r (unstable only)

    tol = 0.00001;

    figure();
    set(gcf,'Position',[100 100 560*1.4 420]) 
    hold on;
    
    % reference idx in Nvec
    i = 1;%length(results_Nvary(end).N);

    % Retrieve modes stability
    unstable_idx = results_Nvary(i).stability_analysis.unstable_idx;
    % nr plots
    nr_plots = min(config.nEigmod, sum(unstable_idx));
    % initialize ylimR
    ylimR = 0;
    
    % Select first or last eigenmodes
    if strcmp(config.whichEigmod,'last')
        start_idx = max(1,sum(unstable_idx)-config.nEigmod)+1;
        end_idx = sum(unstable_idx);
    elseif strcmp(config.whichEigmod,'first')
        start_idx = 1;
        end_idx = nr_plots;
    end
    % Plot eigenmodes
    for j = start_idx:end_idx
        
        % Extract appropriate component
        N = results_Nvary(i).N;
        if strcmp(component, 'V')
            y_data = results_Nvary(i).V(1:N, j);
            output_file = io.outputFile_fig_eigv;
        else % W component
            y_data = results_Nvary(i).V(N+1:2*N, j);
            output_file = io.outputFile_fig_eigw;
        end
        % compute ylim
        ymax = max(y_data);
        last_above_index  =  find(y_data > tol*ymax, 1,'last' ); % index of last point above threshold
        ylimR = max(ylimR, results_Nvary(i).r(last_above_index + 1));
        
        % Set color
        color = colorVec_default(mod(j-start_idx, size(colorVec_default, 1)) + 1, :);
        
        % Plot with appropriate labeling
        plot(results_Nvary(i).r, y_data, 'LineStyle', '-', 'Color', color, ...
             'DisplayName', sprintf('mode $%d$', j-1), 'LineWidth', 1.5);  %-1 since k starts from 0
    end

    % --- COMPARE V, W WITH STATIONARY STATES OF FULL KGW ----------------
    statCmp(config.nEigmod,1)=struct; 
    for j=start_idx:end_idx
        i = j-start_idx+1;
        n = j-1;
        % TEMP: load stationary states of full KGW for comparison
        temp = load("InputStationary/results_k"+sprintf( '%03d', n )+".mat",'eigval','r','f','phi','nSens');
        statCmp(i).n = j-1;
        statCmp(i).mu2 = temp.eigval{temp.nSens};
        statCmp(i).r_prev = temp.r{temp.nSens};
        statCmp(i).u_prev = temp.f{temp.nSens};
        statCmp(i).phi_prev = temp.phi{temp.nSens};
        % Compute u*r, prev
        statCmp(i).ur_prev = statCmp(i).u_prev .* statCmp(i).r_prev;
        % Rescale radius
        [rmax,fmax,~]=findLocMaxima(statCmp(i).r_prev,abs(statCmp(i).ur_prev),n+1);
        rFirstMaxUr = rmax(n+2);
        fFirstMaxUr = fmax(n+2);
        [rmax,fmax,~]=findLocMaxima(results_Nvary(1).r,abs(results_Nvary(1).V(1:results_Nvary(1).N, j)),n+1); % W
        rFirstMaxW = rmax(n+2);
        fFirstMaxW = fmax(n+2);
        r_prev_scaled  = statCmp(i).r_prev / rFirstMaxUr * rFirstMaxW;
        ur_prev_scaled = statCmp(i).ur_prev / fFirstMaxUr * fFirstMaxW;
        % Compute u*r (ok with generic discretization, here use 1)
        statCmp(i).ur   = interp1(r_prev_scaled, ur_prev_scaled,   results_Nvary(1).r, 'spline');
    end
    % OPT: plot stationary states of full KGW for comparison
    if strcmp(component, 'V') && strcmp(config.cmpEigmod, 'on')
        for k=1:nr_plots
            color = [colorVec_default(mod(k-1, size(colorVec_default, 1)) + 1, :), 0.3];
            plot(results_Nvary(1).r,statCmp(k).ur,'--','Color',color, 'LineWidth', 1.5, 'DisplayName',sprintf('$ru_{%d}$',statCmp(k).n));
        end
    end
    %-------------------------------------------------------------------------------------------------------------------------
   

    % Format plot
    xlim([0,ylimR])
    xlabel('$r$');
    ylabel(sprintf('$%s$', component));
    grid on;
    leg = legend('Location', 'northeastoutside'); %,'NumColumns',2);
    set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
    
    if saveFig
        saveas(gcf, io.outputFolder_fig + "/" + output_file);
        fprintf('  Saved: %s\n', output_file);
    end
end

function printStabilityAnalysis(results_Nvary,config)
    % Print summary of stability analysis to command window
    fprintf(sprintf("\n=== Stability Analysis, n=%d ===\n",config.n))
    for i=1:length(results_Nvary)
        fprintf('  N=%d: %d modes, %d stable, %d unstable, %d zero\n', ...
        results_Nvary(i).N, results_Nvary(i).stability_analysis.nr_tot_modes, ...
        results_Nvary(i).stability_analysis.nr_stable_modes, results_Nvary(i).stability_analysis.nr_unstable_modes, results_Nvary(i).stability_analysis.nr_zero_modes);
    end
end


function createSubfolder(folderPath)
    % Create subfolder if it doesn't exist
    if ~isfolder(folderPath)
        mkdir(folderPath);
    end
end