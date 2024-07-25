function [Plot] = plottingMethods(params, FEM, Modal)
    fprintf('%s Plotting methods\n', repmat('-', 1, 4));
    Plot = getPlottingInstance(params, Modal);

    % Plotting the reconstructed FE mesh of the model
    Plot.finiteElementsModel(FEM);

    % For generating 3D animation of vibrations pass "true" as the second argument
    % e.g., Plot.modalExapnsion(Modal, true)
    Plot.modalExpansion(Modal, true)
    Plot.dampedSteadyStateResponse(Modal);
    Plot.dampedTransientResponse(Modal);
    Plot.dampedForcedResponse(Modal);
    Plot.dampedFreeResponse(Modal);

    Plot.impedanceMatrixResponse(Impedance);
    Plot.compareModalAndImpedanceResults(Modal, Impedance);
end


function [Plot] = getPlottingInstance(params, Modal)
    % Create an instance of Plotting class only if it doesn't already exist
    fprintf('%s Creating an instance of the Plotting class\n', repmat('-', 1, 6));
    persistent plot;
    if isempty(plot)
        plot = Plotting(params.plotting_modes, Modal);
    else
        fprintf('%s Already created - skip\n', repmat('-', 1, 8))
    end
    Plot = plot;
end


