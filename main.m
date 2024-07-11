%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main.m
% Main execution file for vibrational analysis of finite elements model of 
% structures
% Initializes variables, manages dependencies, and executes main functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Adjust Default Settings
% Clear workspace, command window, and close all figures
clear; clc; close all; 

t_start = tic;
% Variable format and warning settings
format long;
warning('off', 'MATLAB:MKDIR:DirectoryExists');
% Set default interpreter to LaTeX for all text, axes tick labels, and legends
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');


%% Initialization

disp('* Start vibrational analysis of finite elements model of structures...');
% Add project folder containig files for model parameters and load them
addpath('Inputs', "Library");
params = load('parameters.mat');
% Add number of low frequency modes to consider for plotting mode shapes
params.low_modes = 6;
% Number of excitation frequencies (w) to be interpolated in the range of significant modes
params.num_w = 150;
% Modes to plot
params.plotting_modes = [10, 15];

% Initializing a parallel pool (used for computing impedance matrix response)
if isempty(gcp('nocreate'))
    parpool;
end

% Initializing class instances for the model and each method
disp('-- Initializing...');
% Finite elements model (FEM)
FEM = getObjectInstances("FEM", params); 
% Modal expansion method (Modal)
Modal = getObjectInstances("Modal", params);
% Impedance matrix method (Impedance) 
Impedance = getObjectInstances("Impedance", params);
% Plotting (Plot) 
Plot = getObjectInstances("Plot", params);


%% Modeling and Analysis

% Modal expansion method
% Compute the damped steady-state and transient responses for harmonic excitations
Modal = Modal.computeForcedHarmonicResponse(params, FEM);
% Compute the damped free response with final value of forced response as its initial condition
Modal = Modal.computeDampedFreeResponse(params);

% Impedance matrix method
% Compute the damped steady-state response using the impedance matrix method
Impedance = Impedance.computeForcedHarmonicResponse(params, FEM);


%% Plotting

% For generating 3D animation of vibrations pass "true" as the second argument 
% e.g., Plot.modalExapnsion(Modal, true)
Plot.modalExpansion(Modal)
Plot.dampedSteadyStateResponse(Modal);
Plot.dampedTransientResponse(Modal);
Plot.dampedForcedResponse(Modal);
Plot.dampedFreeResponse(Modal);
Plot.impedanceMatrixResponse(Impedance);
Plot.compareModalAndImpedanceResults(Modal, Impedance);

% Calculate and display the execution time
delta_t = seconds(toc(t_start));
delta_t.Format = 'hh:mm:ss.SSS';
fprintf('* Successfully completed the execution in %s *\n', char(delta_t));
