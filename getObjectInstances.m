function [FEM, Modal, Impedance, Plot] = getObjectInstances(params)
    % Create the objects instances only if they don't already exist
    disp('-- Initializing...');
    persistent fem modal impedance plot;
    
    if isempty(fem)
        fem = FiniteElementsModel(params);
    else
        disp('---- Finite elements model already initialized - skip')
    end
    FEM = fem;
    
    if isempty(modal)
        modal = ModalExpansionMethod(params, FEM);
    else
        disp('---- Modal expansion method already initialized - skip')
    end
    Modal = modal;
    
    if isempty(impedance)
        impedance = ImpedanceMatrixMethod(params, FEM);
    else
        disp('---- Impedance matrix method already initialized - skip')
    end
    Impedance = impedance;
    
    if isempty(plot)
        plot = Plotting();
    else
        disp('---- Plotting already initialized - skip')
    end
    Plot = plot;
end