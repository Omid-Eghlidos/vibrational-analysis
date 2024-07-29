%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main.m
% Main execution file for vibrational analysis of FE model of structures
% Initializes variables, manages dependencies, and executes main functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
disp('* Start vibrational analysis of finite elements model of structures...');

%% Initialization
% Adjust the results to be plotted (1 = true, 0 = false)
plotting.FEM       = struct('model',        1);
plotting.Modal     = struct('frequencies',  1, 'low_modes',  1, ...
                            'animation',    1, 'Usteady',    1, ...
                            'Utransient',   1, 'Uforced',    1, ...
                            'Ufree',        1, 'modes2plot', [10, 15]);
plotting.Impedance = struct('Usteady',      1, 'modes2plot', [10, 15], ...
                            'versus_modal', 1);
plotting.CB        = struct('frequencies',  1);
% Initialize the rest of parameters
params = Initialize(plotting);


%% Analysis
% Set the type of analysis to run (1 = true, 0 = false)
run = struct('Modal', 1, 'Impedance', 1, 'Craig_Bampton', 1, 'Plotting', 1);
results = runAnalysis(params, run);
