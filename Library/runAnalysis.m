function [results] = runAnalysis(params, run)
    % Will run the specified analaysis to obtain results
    t_start = tic;
    results = struct();
    disp('-- Starting analysis')
    
    FEM = finiteElementsMethod(params.FEM);
    if run.Plotting
        plottingMethods(params.plotting, FEM);
    end
    results.FEM = FEM;
    
    if run.Modal
        Modal = modalExpansionMethod(params.Modal, FEM);
        if run.Plotting
            plottingMethods(params.plotting, Modal);
        end
        results.Modal = Modal;
    end

    % Impedance only runs if Modal has already executed
    if run.Impedance && isfield(results, 'Modal')
        Impedance = impedanceMatrixMethod(params.Impedance, FEM, Modal);
        if run.Plotting
            plottingMethods(params.plotting, Impedance);
            plottingMethods(params.plotting, Impedance, Modal)
        end
        results.Impedance = Impedance;
    end

    if run.Craig_Bampton
        CB = craigBamptonMethod(params.CB, FEM);
        if run.Plotting
            plottingMethods(params.plotting, CB);
        end
        results.CB = CB;
    end
    
    % Runtime statistics
    delta_t = seconds(toc(t_start));
    delta_t.Format = 'hh:mm:ss.SSS';
    fprintf('* Execution completed in %s *\n', char(delta_t));
    clear run;
end


function [FEM] = finiteElementsMethod(params)
    % Read the FEM model matrices from files and create the structure's mesh 
    fprintf('%s Finite elements method\n', repmat('-', 1, 4));
    FEM = getFiniteElementsModelInstance(params);
end


function [Modal] = modalExpansionMethod(params, FEM)
    % Compute the damped forced and free response using modal expansion method
    fprintf('%s Modal expansion method\n', repmat('-', 1, 4));
    Modal = getModalExpansionInstance(params, FEM);
    Modal = Modal.computeForcedHarmonicResponse(params, FEM);
    Modal = Modal.computeDampedFreeResponse(params);
end


function [Impedance] = impedanceMatrixMethod(params, FEM, Modal)
    % Compute the damped forced response using impedance matrix method
    fprintf('%s Impedance matrix method\n', repmat('-', 1, 4));
    Impedance = getImpedanceMatrixInstance(Modal);
    Impedance = Impedance.computeForcedHarmonicResponse(params, FEM);
end


function [CB] = craigBamptonMethod(params, FEM)
    % Perform Craig-Bampton reduction by selecting substructures on the mesh
    fprintf('%s Craig-Bampton method\n', repmat('-', 1, 4));
    CB = getCraigBamptonInstance(params);
    CB = CB.addSubstructures(FEM);
    CB = CB.performCraigBamptonReduction();
end


function plottingMethods(params, varargin)
    % Depending on the provided input plot the results of the analysis
    Plot = getPlottingInstance();
    % Determine the type of input
    if length(varargin) == 1
        arg = varargin{1};
        if class(arg) == "FiniteElementsModel"
            Plot.finiteElementsModel(params.FEM, arg);
        elseif class(arg) == "ModalExpansion"
            Plot.modalExpansionMethod(params.Modal, arg);
        elseif class(arg) == "ImpedanceMatrix"
            Plot.impedanceMatrixMethod(params.Impedance, arg);
        elseif class(arg) == "CraigBampton"
            Plot.craigBamptonMethod(params.CB, arg);
        end
    elseif length(varargin) == 2
        Plot.impedanceMatrixMethod(params.Impedance, varargin{1}, varargin{2});
    end
end


function [FEM] = getFiniteElementsModelInstance(params)
    % Create an instance of FiniteElementsModel class only if it doesn't already exist
    persistent fem;
    if isempty(fem)
        fem = FiniteElementsModel(params);
    end
    FEM = fem;
end


function [Modal] = getModalExpansionInstance(params, FEM)
    % Create an instance of ModalExpansion class only if it doesn't already exist
    persistent modal;
    if isempty(modal)
        modal = ModalExpansion(params, FEM);
    end
    Modal = modal;
end


function [Impedance] = getImpedanceMatrixInstance(Modal)
    % Create an instance of ImpedanceMatrix class only if it doesn't already exist
    persistent impedance;
    if isempty(impedance)
        impedance = ImpedanceMatrix(Modal);
    end
    Impedance = impedance;
end


function [Craig_Bampton] = getCraigBamptonInstance(params)
    % Create an instance of CraigBampton class only if it doesn't already exist
    persistent craig_bampton;
    if isempty(craig_bampton)
        craig_bampton = CraigBampton(params);
    end
    Craig_Bampton = craig_bampton;
end


function [Plot] = getPlottingInstance()
    % Create an instance of Plotting class only if it doesn't already exist
    persistent plot;
    if isempty(plot)
        plot = Plotting();
    end
    Plot = plot;
end
