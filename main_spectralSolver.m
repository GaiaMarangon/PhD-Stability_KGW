%----------------------------------------------------------------------
% Main script for solving KGW spectral problems with convergence analysis
%----------------------------------------------------------------------
clear all, close all, clc

% Physical Parameters
params.n = 3;                  % Excitation index for background stationary state

% Numerical parameters
% params.Nvec = [200,250,300,350,400,450,500,550,600];    % Number of discretization points
params.Nvec = linspace(1000,1500,11)';    % Number of discretization points
params.tolerance = 1e-10;      % Tolerance for zero modes
params.max_eigenvalue_cutoff = []; % Threshold for spurious modes (optional)
params.max_heuristic_cutoff = []; % Percentage of discarded large eigenvalues (optional, heuristic)

% Input/Output parameters
io.inputFolder = "InputStationary";
io.outputFolder_res = "Results";
% Results files
io.outputFile_res_base = "KGWspectral" + sprintf('_n%03d', params.n);
io.outputFile_res_all = io.outputFile_res_base + ".mat";

% Create output folders
createSubfolder("./" + io.outputFolder_res);

% Load background data
[mu2_n, r_prev, u_prev, phi_prev] = loadBackgroundData(io.inputFolder, params.n);
params.r0 = r_prev(1);
params.rmax = r_prev(end);

% Initialize or load existing compatible results
results_file = io.outputFolder_res + "/" + io.outputFile_res_all;
[results_Nvary, existing_params, is_new_dataset] = initializeOrLoadResults(results_file, params);

% Main computation loop - iterate over discretization points
fprintf('\n=== STARTING KGW SPECTRAL ANALYSIS ===\n');
total_calculations = 0;
for i = 1:length(params.Nvec)
    params.currentN = params.Nvec(i);
    
    % Check if this specific N already exists in the dataset
    existing_idx = findExistingN(results_Nvary, params.currentN);
    
    if existing_idx > 0
        fprintf('N = %d (%d/%d): Already computed, skipping\n', params.currentN, i, length(params.Nvec));
        continue;
    end
    
    fprintf('\nSolving for N = %d (%d/%d)\n', params.currentN, i, length(params.Nvec));
    total_calculations = total_calculations + 1;
    
    % Solve spectral problem
    new_result = solveKGWPerturbations(mu2_n, r_prev, u_prev, phi_prev, params);
    
    % Find correct position to insert (keep N sorted)
    insert_pos = findInsertPosition(results_Nvary, params.currentN);
    results_Nvary = insertResult(results_Nvary, new_result, insert_pos);
    
    % Print summary
    printSpectralSummary(new_result.stability_analysis, params.currentN);
    
    % Save intermediate results after each calculation    
    save(results_file, 'results_Nvary', 'params', 'mu2_n');
end

% Final save
save(results_file, 'results_Nvary', 'params', 'mu2_n');

% Print summary
printComputationSummary(results_file,results_Nvary,params,total_calculations)



%----------------------------------------------------------------------
%--- SUPPORT FUNCTIONS -----------------------------------------------
%----------------------------------------------------------------------

function [mu2_n, r_prev, u_prev, phi_prev] = loadBackgroundData(inputFolder, n)
    % Load background stationary states from input folder
    filename = inputFolder + "/results_k" + sprintf('%03d', n) + ".mat";
    if ~isfile(filename)
        error('Background data file not found: %s', filename);
    end
    
    load(filename, 'eigval', 'r', 'f', 'phi', 'nSens');
    mu2_n = eigval{nSens};
    r_prev = r{nSens};
    u_prev = f{nSens};
    phi_prev = phi{nSens};
    
    fprintf('Loaded background data: n=%d\n\n', n);
end

function [results_Nvary, existing_params, is_new_dataset] = initializeOrLoadResults(results_file, params)
    % Initialize results array or load existing compatible dataset
    
    % Create new dataset
    if ~isfile(results_file)
        fprintf('No existing results found. Starting new dataset.\n');
        results_Nvary = [];
        existing_params = params;
        is_new_dataset = true;
        return;
    end
    
    % Load existing dataset
    existing_data = load(results_file);
    
    % Check if existing data is compatible (same parameters except Nvec)
    if isDatasetCompatible(existing_data, params)
        fprintf('Found compatible existing dataset. Will add missing N values.\n');
        results_Nvary = existing_data.results_Nvary;
        existing_params = existing_data.params;
        is_new_dataset = false;
        
        % Merge and sort all N values (existing + requested)
        stored_Nvec = existing_data.params.Nvec;
        unified_Nvec = unique([stored_Nvec, params.Nvec]);
        existing_params.Nvec = unified_Nvec;
        
    else
        fprintf('Existing dataset has incompatible parameters. Starting new dataset. Will overwrite old parameters.\n\n');
        disp('Do you wish to continue? Press a key to continue, Ctrl+C to quit.')  
        pause;
        results_Nvary = [];
        existing_params = params;
        is_new_dataset = true;
    end
end

function is_compatible = isDatasetCompatible(existing_data, current_params)
    % Check if existing dataset is compatible with current parameters
    % Compatible means all parameters are the same except Nvec
    
    is_compatible = false;
    
    % Check critical parameters (all except Nvec)
    fields_to_check = {'n', 'tolerance', 'r0', 'rmax','max_eigenvalue_cutoff', 'max_heuristic_cutoff'};
    for i = 1:length(fields_to_check)
        field = fields_to_check{i};
        if ~isfield(existing_data.params, field) || ...
           ~isequal(existing_data.params.(field), current_params.(field))
            fprintf('  Parameter mismatch in %s\n', field);
            fprintf('- Old value: %s \n- New value: %s \n', [num2str(existing_data.params.(field)), num2str(current_params.(field))])
            return;
        end
    end
    
    is_compatible = true;
end

function existing_idx = findExistingN(results_Nvary, target_N)
    % Find if target_N already exists in the results array
    % Returns index if found, 0 if not found
    
    existing_idx = 0;
    if isempty(results_Nvary)
        return;
    end
    
    for i = 1:length(results_Nvary)
        if results_Nvary(i).N == target_N
            existing_idx = i;
            return;
        end
    end
end

function insert_pos = findInsertPosition(results_Nvary, new_N)
    % Find the correct position to insert new result to keep N values sorted
    
    if isempty(results_Nvary)
        insert_pos = 1;
        return;
    end
    
    % Find position where new_N should be inserted
    for i = 1:length(results_Nvary)
        if new_N < results_Nvary(i).N
            insert_pos = i;
            return;
        end
    end
    
    % If we get here, new_N is larger than all existing N values
    insert_pos = length(results_Nvary) + 1;
end

function results_Nvary = insertResult(results_Nvary, new_result, insert_pos)
    % Insert new result at the specified position, keeping array sorted by N
    
    if isempty(results_Nvary)
        results_Nvary = new_result;
    else
        % Shift elements and insert
        results_Nvary = [results_Nvary(1:insert_pos-1), new_result, results_Nvary(insert_pos:end)];
    end
end

function [new_result] = solveKGWPerturbations(mu2_n, r_prev, u_prev, phi_prev, params)
    % Solve radial KGW linearized perturbation problem
    
    % Extract parameters
    N = params.currentN;
    rMax = params.rmax;
    r0 = params.r0;
    
    % Define uniform grid and differentiation matrices
    [r, D2, u_n, phi_n] = createUniformFiniteDiffGrid(N, r0, rMax, r_prev, u_prev, phi_prev);
    
    % Construct and solve spectral matrix
    A = constructSpectralMatrix(mu2_n, u_n, phi_n, D2);
    A = applyBoundaryConditions(A);
    
    % Solve eigenvalue problem
    [V, E] = eig(A);
    omega2 = diag(E);
    
    % Sort by eigenvalue
    [omega2, idx] = sort(omega2);
    V = V(:, idx);
    
    % Filter spurious modes if requested
    [omega2, V] = filterSpuriousModes(omega2, V, params);
    
    % Normalize eigenvectors and add boundary conditions
    V = normalizeEigenvectors(V, r);
    
    % Analyze stability
    stability_analysis = analyzeStability(omega2, params);

    % Store results
    new_result = struct();
    new_result.N = params.currentN;
    new_result.omega2 = omega2;
    new_result.V = V;
    new_result.stability_analysis = stability_analysis;
    new_result.r = r;
    new_result.u_n = u_n;
    new_result.phi_n = phi_n;
end

function [r, D2, u_n, phi_n] = createUniformFiniteDiffGrid(N, r0, rMax, r_prev, u_prev, phi_prev)
    % Create uniform grid with finite difference second derivative matrix
    
    % Uniform grid points
    r = linspace(r0, rMax, N)';
    
    % Second derivative matrix - centered finite differences
    dr = (rMax - r0) / (N - 1);
    D2 = (diag(ones(N-1,1), 1) - 2*eye(N) + diag(ones(N-1,1), -1)) / dr^2;
    
    % Interpolate background functions
    u_n = interp1(r_prev, u_prev, r, 'spline');
    phi_n = interp1(r_prev, phi_prev, r, 'spline');
    
    % Verify matrix properties
    sym_error = max(max(abs(D2 - D2')));
    if sym_error > 1e-14
        warning('D2 matrix symmetry error: %e', sym_error);
    end
end

function A = constructSpectralMatrix(mu2_n, u_n, phi_n, D2)
    % Construct the spectral matrix for KGW perturbation problem
    % | mu^2 + 2*phi_n - d^2/dr^2    2*u_n     | |V|     |V|
    % |        2*u_n              -d^2/dr^2    | |W| = ω²|W|
    
    A11 = diag(mu2_n + 2*phi_n) - D2;
    A12 = 2*diag(u_n);
    A22 = -D2;
    
    A = [A11, A12;
         A12, A22];
end

function A = applyBoundaryConditions(A)
    % Apply Dirichlet boundary conditions
    N = size(A, 1) / 2;
    
    % Remove first and last rows/columns of each block (Dirichlet BC)
    A = [A(2:N-1, 2:N-1),         A(2:N-1, N+2:2*N-1);
         A(N+2:2*N-1, 2:N-1),     A(N+2:2*N-1, N+2:2*N-1)];
end

function [omega2_clean, V_clean] = filterSpuriousModes(omega2, V, params)
    % Filter spurious modes based on eigenvalue magnitude
    
    if ~isempty(params.max_eigenvalue_cutoff)
        valid_idx = abs(omega2) < params.max_eigenvalue_cutoff;
    elseif ~isempty(params.max_heuristic_cutoff)
        sorted_abs = sort(abs(omega2), 'descend');
        cutoff_idx = max(1, round(params.max_heuristic_cutoff * length(omega2)));
        threshold = sorted_abs(cutoff_idx);
        valid_idx = abs(omega2) <= threshold;
    else
        valid_idx = true(size(omega2));
    end
    
    omega2_clean = omega2(valid_idx);
    V_clean = V(:, valid_idx);
    
    n_removed = sum(~valid_idx);
    if n_removed > 0
        fprintf('  Removed %d spurious modes\n', n_removed);
    end
end

function V_normalized = normalizeEigenvectors(V, r)
    % Normalize eigenvectors with L2 norm and add boundary conditions
    
    N = size(V, 1) / 2 + 2;
    V_normalized = zeros(2*N, size(V, 2));
    
    for i = 1:size(V, 2)
        % Extract components
        v = V(1:(N-2), i);
        w = V((N-2)+1:end, i);
        
        % Sign convention
        if v(1) < 0, v = -v; end
        if w(1) < 0, w = -w; end
        
        % L2 normalization with grid spacing
        dr = r(2:end-1) - r(1:end-2);
        norm_v = sqrt(sum(dr .* (v.^2)));
        norm_w = sqrt(sum(dr .* (w.^2)));
        
        % Store normalized vectors with boundary conditions
        V_normalized(2:N-1, i) = v / norm_v;
        V_normalized(N+2:end-1, i) = w / norm_w;
    end
end

function stability = analyzeStability(omega2, params)
    % Analyze stability by counting stable/unstable modes
    
    % Indices of stable, unstable and zero modes
    stability.stable_idx = omega2 >= params.tolerance;
    stability.unstable_idx = omega2 <= -params.tolerance;
    stability.zero_idx = abs(omega2) < params.tolerance;
    % Total number of stable, unstable and zero modes
    stability.nr_tot_modes = length(omega2);
    stability.nr_stable_modes = sum(stability.stable_idx);
    stability.nr_unstable_modes = sum(stability.unstable_idx);
    stability.nr_zero_modes = sum(stability.zero_idx);
end

function printSpectralSummary(stability, N)
    % Print summary of spectral analysis
    fprintf('  N=%d: %d modes, %d stable, %d unstable, %d zero\n', ...
        N, stability.nr_tot_modes, stability.nr_stable_modes, stability.nr_unstable_modes, stability.nr_zero_modes);
end

function createSubfolder(folderPath)
    % Create subfolder if it doesn't exist
    if ~isfolder(folderPath)
        mkdir(folderPath);
    end
end

function printComputationSummary(results_file,results_Nvary,params,total_calculations)
    % Print summary of the computation

    % Extract all N values from results_Nvary struct 
    computed_N_values = arrayfun(@(x) x.N, results_Nvary);
    % Print
    fprintf('\n------------------------------------------------------------\n')
    fprintf('\n=== COMPUTATION COMPLETED ===\n');
    
    fprintf('\nDataset n=%d contains N values: [%s]\n', [params.n, sprintf('%d ', computed_N_values)]);
    if total_calculations == 0
        fprintf('All requested calculations were already present in the dataset.\n');
    else
        fprintf('Added %d new calculations to the dataset.\n', total_calculations);
    end
    
    fprintf('\nSaving final results to: %s\n', results_file);
end