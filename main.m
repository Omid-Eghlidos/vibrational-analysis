%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main.m
% Main execution file for vibrational analysis of finite elements model of 
% structures
% Initializes variables, manages dependencies, and executes main functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Adjust Default Settings
% Clear workspace, command window, and close all figures
clear; clc; close all; 

disp('* Start vibrational analysis of finite elements model of structures...');
t_start = tic;
% Variable format and warning settings
format long;
warning('off', 'MATLAB:MKDIR:DirectoryExists');
% Set default interpreter to LaTeX for all text, axes tick labels, and legends
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');


%% Initialization

disp('-- Initializing...');
% Add project folder containig files for model parameters and load them
addpath('Inputs', 'Library');
params = load('parameters.mat');
% Add number of low frequency modes to consider for plotting mode shapes
params.low_modes = 6;
% Number of excitation frequencies (w) to be interpolated in the range of significant modes
params.num_w = 150;
% Modes to plot
params.plotting_modes = [10, 15];

%% Modeling and Analysis

disp('-- Starting analysis')
% Finite elements model (FEM)
FEM = finiteElementsMethod(params);
% Modal expansion method
%Modal = modalExpansionMethod(params, FEM);
% Impedance matrix method
%Impedance = impedanceMatrixMethod(params, FEM, Modal);
% Craig-Bampton method
CB = craigBamptonMethod(FEM);
% Plotting
%plottingMethods(params, FEM, Modal)


%% Running statistics
% Calculate and display the execution time
delta_t = seconds(toc(t_start));
delta_t.Format = 'hh:mm:ss.SSS';
fprintf('* Execution completed in %s *\n', char(delta_t));
