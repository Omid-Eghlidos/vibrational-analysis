%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main.m
% Main execution file for vibrational analysis of finite elements model of 
% structures
% Initializes variables, manages dependencies, and executes main functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
% Clear workspace, command window, and close all figures
clear; clc; close all; 

% Variable format and warning settings
format long;
warning('off', 'MATLAB:MKDIR:DirectoryExists');
% Set default interpreter to LaTeX for all text, axes tick labels, and legends
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

disp('* Start vibrational analysis of a finite element model of structures...');

% Add project folder containig files for model parameters and load them
addpath('inputs')
params = load('project26.mat');
% Add number of low frequency modes to consider for plotting mode shapes
params.low_modes = 6;

% Initializing finite elements model (FEM), modal expansion method (Modal)
% impedance matrix method, and plotting parameters
[FEM, Modal, Impedance, Plot] = getObjectInstances(params);


%% Modeling and Analysis

% Modal expansion method
% Compute the damped steady-state and transient responses due to harmonic excitations
%Modal = Modal.computeForcedHarmonicResponse(params, FEM);
% Compute the damped free response of the system
%Modal = Modal.computeDampedFreeResponse(params);

% Impedance matrix method
% Compute the damped steady-state response using Impedance matrix
Impedance.computeForcedHarmonicResponse(params, FEM);


%% Plotting

% Plot.modalExpansion(Modal)
% Plot.dampedSteadyStateResponse(Modal);
% Plot.dampedTransientResponse(Modal);
% Plot.dampedForcedResponse(Modal);
% Plot.dampedFreeResponse(Modal);
Plot.impedanceMatrixResponse(Impedance);

disp('* Successful execution *');
