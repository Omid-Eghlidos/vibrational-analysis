function [params] = Initialize(plotting)
    % Initialize required parameters and settings for running analyses
    % Load parameters and settings if it has not already
    disp('-- Initializing');
    % Variable format and warning settings
    format long;
    warning('off', 'MATLAB:MKDIR:DirectoryExists');
    % Set default interpreter to LaTeX for all text, axes tick labels, and legends
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');    
    % Add input folder containig files for model parameters
    addpath('Inputs', 'Library');

    % Finite elements model parameters
    % Rayleigh damping coefficients of alpha and beta
    params.FEM.alpha = 30.2478002081387;
    params.FEM.beta  = 7.748497358063650e-07;
    % Harmonic force: [node, DoF, magnitude]
    params.FEM.Fc = [94, 2, 1];
    params.FEM.Fs = [];

    % Modal expansion method parameters
    % Rayleigh damping coefficients of alpha and beta
    params.Modal.alpha = 30.2478002081387;
    params.Modal.beta  = 7.748497358063650e-07;
    % Harmonic cos and sin forces: [node, DoF, magnitude]
    params.Modal.Fc = [94, 2, 1];
    params.Modal.Fs = [];
    % Range of frequencies within which computing harmonic response
    params.Modal.range_harmonic = [0; 1000];
    % Range of frequencies within which computing free response
    params.Modal.free_range = [0; 1000];
    % Nodes and their corresponding DoF to compute the harmonic response 
    params.Modal.nodes_harmonic = [94, 2; 204, 2; 166, 2; 238, 2];
    % Nodes and their corresponding DoF to compute the transient response
    params.Modal.nodes_transient = [94, 2; 204, 2; 166, 2; 238, 2];
    % Range of frequencies and their corresponding damping ratios be interpolated
    params.Modal.zetas = [0, 0.01; 500, 0.01; 1000, 0.05; 100000000, 0.05];
    % Time (s) up to which computing the steady-state and transient response
    params.Modal.time0 = 0.30;
    % Time (s) after time0 up to which compute the free damped response
    params.Modal.time1 = 0.60;
    % Number of low frequency modes for plotting and making 3D animations
    params.Modal.low_modes = 6;
    % Number of excitation frequencies (w) to be interpolated 
    % NOTE: Higher number increases runtime specially for Impedance method
    params.Modal.num_w = 150;

    % Impedance matrix method parameters
    % Harmonic force: [node, DoF, magnitude]
    params.Impedance.Fc = [94, 2, 1];
    params.Impedance.Fs = [];
    % Nodes and their corresponding DoF to compute the harmonic response 
    params.Impedance.nodes_harmonic = [94, 2; 204, 2; 166, 2; 238, 2];
    
    % Craig-Bampton method settings
    params.CB.use_existing_substructures = 0;
    params.CB.tolerance = 1e-3;
    params.CB.output_path = "Results/Craig_Bampton";
    params.CB.Nmodes = 15;
    
    % Plotting settings
    params.plotting = plotting;
    clear plotting;
end
