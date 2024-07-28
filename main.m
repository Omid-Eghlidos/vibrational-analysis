%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main.m
% Main execution file for vibrational analysis of FE model of structures
% Initializes variables, manages dependencies, and executes main functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
disp('* Start vibrational analysis of finite elements model of structures...');

% Plotting settings (1 = true, 0 = false)
plotting.FEM = struct('model', 1);
plotting.Modal = struct('frequencies', 1, 'low_modes', 1, 'animation', 0, ...
    'modes2plot', [10, 15],'Usteady', 1, 'Utransient', 1, 'Uforced', 1, 'Ufree', 1);
plotting.Impedance = struct('Usteady', 1, 'versus_modal', 1);
plotting.CB = struct('frequencies', 1);
% Analysis to run
run = struct('Modal', 1, 'Impedance', 1, 'Craig_Bampton', 1, 'Plotting', 1);

% Initialization
params = Initialize(plotting);
% Modeling and Analysis
analysis = runAnalysis(params, run);


function [params] = Initialize(plotting)
    % Initialize required parameters and settings
    % Load parameters and settings if it has not already
    if exist('params', 'var') ~= 0
        return;
    end
    disp('-- Initializing');
    % Variable format and warning settings
    format long;
    warning('off', 'MATLAB:MKDIR:DirectoryExists');
    % Set default interpreter to LaTeX for all text, axes tick labels, and legends
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    
    % Add project folder containig files for model parameters and load them
    addpath('Inputs', 'Library');
    params = load('parameters.mat');
    % Modal expansion method settings
    % Number of low frequency modes for plotting
    params.low_modes = 6;
    % Number of excitation frequencies (w) to be interpolated
    params.num_w = 100;
    
    % Craig-Bampton method settings
    params.CB.use_existing_substructures = 1;
    params.CB.tolerance = 1e-3;
    params.CB.output_file = "Inputs/substructures.mat";
    params.CB.Nmodes = 15;
    
    % Plotting settings
    params.plotting = plotting;
    clear plotting;
end


function [analysis] = runAnalysis(params, run)
    % Depending on the selection will run any of the analaysis
    t_start = tic;
    analysis = struct();
    disp('-- Starting analysis')
    persistent FEM Modal Impedance CB;
    
    FEM = finiteElementsMethod(params);
    if run.Plotting
        plottingMethods(params.plotting, FEM);
    end
    analysis.FEM = FEM;
    
    if run.Modal
        Modal = modalExpansionMethod(params, FEM);
        if run.Plotting
            plottingMethods(params.plotting, Modal);
        end
        analysis.Modal = Modal;

        % Impedance only runs if Modal has already executed
        if run.Impedance
            Impedance = impedanceMatrixMethod(params, FEM, Modal);
            if run.Plotting
                plottingMethods(params.plotting, Impedance);
                plottingMethods(params.plotting, Impedance, Modal)
            end
            analysis.Impedance = Impedance;
        end
    end

    if run.Craig_Bampton
        CB = craigBamptonMethod(params.CB, FEM);
        if run.Plotting
            plottingMethods(params.plotting, CB);
        end
        analysis.CB = CB;
    end
    
    % Runtime
    delta_t = seconds(toc(t_start));
    delta_t.Format = 'hh:mm:ss.SSS';
    fprintf('* Execution completed in %s *\n', char(delta_t));
    clear run;
end
