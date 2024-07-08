%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main.m
% Main execution file for vibrational analysis of finite elements model of 
% structures
% Initializes variables, manages dependencies, and executes main functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Adjust Default Settings
% Clear workspace, command window, and close all figures
clear; clc; close all; 

% Variable format and warning settings
format long;
warning('off', 'MATLAB:MKDIR:DirectoryExists');
% Set default interpreter to LaTeX for all text, axes tick labels, and legends
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

%% Initialization
disp('* Start vibrational analysis of finite element model of a structure...');

% Add project folder containig files for model parameters and load them
addpath('inputs')
params = load('parameters.mat');
% Add number of low frequency modes to consider for plotting mode shapes
params.low_modes = 6;
% Initializing finite elements model (FEM), modal expansion method (Modal)
% impedance matrix method, and plotting parameters
[FEM, Modal, Impedance, Plot] = getObjectInstances(params);


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

Plot.modalExpansion(Modal)
Plot.dampedSteadyStateResponse(Modal);
Plot.dampedTransientResponse(Modal);
Plot.dampedForcedResponse(Modal);
Plot.dampedFreeResponse(Modal);
Plot.impedanceMatrixResponse(Impedance);

disp('* Successful execution *');
