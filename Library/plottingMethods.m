function plottingMethods(params, varargin)
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


function [Plot] = getPlottingInstance()
    % Create an instance of Plotting class only if it doesn't already exist
    persistent plot;
    if isempty(plot)
        plot = Plotting();
    end
    Plot = plot;
end
