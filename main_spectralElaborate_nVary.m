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
config.saveFigures = true;      % Flag to control figure saving
config.nVec = linspace(20,30,6);     % Vector of n values to be used for comparison across different n

% Input/Output parameters
io.inputFolder_res = "Results";
io.outputFolder_fig = "Figures";
io.inputFile_res = @(n) "KGWspectral" + sprintf('_n%03d', n) + ".mat";

% Figure output files
io.outputFile_fig_specIdx_nVary = "KGWspectral_nVary_spectrum_Idx.png";

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



%--- ANALYSIS ACROSS DIFFERENT n VALUES ---------------------------------

% 1. Plot spectrum as n varies
plotSpectrum_nVary(config, io, config.saveFigures);

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

%----------------------------------------------------------------------
function plotSpectrum_nVary(config, io, saveFig)
    % Plot Spectrum against index as n varies

    figure()
    set(gcf,'Position',[100 100 560*1 420])
    hold on

    for i=1:length(config.nVec)
        
        % Set n
        n = config.nVec(i);

        % Load data 
        ion = io;
        ion.inputFile_res = io.inputFile_res(n);
        [~, results_Nvary] = loadResults(ion);

        % Plot Settings
        alpha = 1 - (i-1)/ (length(config.nVec));
        
        % Reference eigenvalues
        ref_omega2 = results_Nvary(end).omega2 ;
            % % Retrieve stable indices, sort ascending
            % [asc_omega2,idx_asc] = sort(ref_omega2);
            % unstable_idx = results_Nvary(end).stability_analysis.unstable_idx(idx_asc);
            % unstable_idx_progressive = cumsum(unstable_idx);
        % Retrieve stable indices, sort descending
        [desc_omega2,idx_desc] = sort(ref_omega2,'descend');
        unstable_idx = results_Nvary(end).stability_analysis.unstable_idx(idx_desc);
        unstable_idx_progressive = cumsum(unstable_idx);

        % Set xdata, ydata
        xdata_unstab = unstable_idx_progressive(unstable_idx) -1; %-1 since k starts from 0
        ydata_unstab = desc_omega2(unstable_idx);
            % % Power law fits ( log(y) = p(1)*log(x) + p(2) )
            % p_unstab = polyfit( log(xdata_unstab), log(ydata_unstab),1);
    
        % Add plots
            % plot(xdata_unstab, exp(p_unstab(2))* xdata_unstab .^p_unstab(1), '-', 'Color', [colorVec_default(1,:) alpha], 'DisplayName', sprintf('$\omega^2_{%d} = %+.2f k^%+.2f$',[n,exp(p_unstab(2)),p_unstab(1)]))
        scatter(xdata_unstab, ydata_unstab, 15, 'filled','MarkerFaceColor',config.colorVec_default(1,:),'MarkerFaceAlpha', alpha, 'DisplayName', sprintf('n= %d',n));

    end
    
    % Format Plot
    xlabel('$k$');
    ylabel('$\omega^2$');
    leg = legend('Location', 'southwest');
    set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
    % set(gca,'xscale','log')
    % set(gca,'yscale','log')

    % Save Figure
    if saveFig
        saveas(gcf, io.outputFolder_fig + "/" + io.outputFile_fig_specIdx_nVary);
        fprintf('  Saved: %s\n', io.outputFile_fig_specIdx_nVary);
    end

end

%----------------------------------------------------------------------
function createSubfolder(folderPath)
    % Create subfolder if it doesn't exist
    if ~isfolder(folderPath)
        mkdir(folderPath);
    end
end