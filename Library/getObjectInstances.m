function [Instance] = getObjectInstances(instance, params)
    % Create the objects instances only if they don't already exist
    persistent FEM Modal Impedance Plot;
    if instance == "FEM"
        FEM = getFiniteElementsModelInstance(params);
        Instance = FEM;
    elseif instance == "Modal"
        Modal = getModalExpansionMethodInstance(params, FEM);
        Instance = Modal;
    elseif instance == "Impedance"
        Impedance = getImpedanceMatrixMethodInstance(Modal);
        Instance = Impedance;
    elseif instance == "Plot"
        Plot = getPlottingInstance(params, Modal);
        Instance = Plot;
    end
end


function [FEM] = getFiniteElementsModelInstance(params)
    persistent fem;
    if isempty(fem)
        fem = FiniteElementsModel(params);
    else
        disp('---- Finite elements model already initialized - skip')
    end
    FEM = fem;
end


function [Modal] = getModalExpansionMethodInstance(params, FEM)
    persistent modal;
    if isempty(modal)
        modal = ModalExpansionMethod(params, FEM);
    else
        disp('---- Modal expansion method already initialized - skip')
    end
    Modal = modal;
end


function [Impedance] = getImpedanceMatrixMethodInstance(Modal)
    persistent impedance;
    if isempty(impedance)
        impedance = ImpedanceMatrixMethod(Modal);
    else
        disp('---- Impedance matrix method already initialized - skip')
    end
    Impedance = impedance;
end


function [Plot] = getPlottingInstance(params, Modal)
    persistent plot;
    if isempty(plot)
        plot = Plotting(params.plotting_modes, Modal);
    else
        disp('---- Plotting already initialized - skip')
    end
    Plot = plot;
end
