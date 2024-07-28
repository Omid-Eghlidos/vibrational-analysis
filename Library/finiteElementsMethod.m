function [FEM] = finiteElementsMethod(params)
    fprintf('%s Finite elements method\n', repmat('-', 1, 4));
    FEM = getFiniteElementsModelInstance(params);
end


function [FEM] = getFiniteElementsModelInstance(params)
    % Create an instance of FiniteElementsModel class only if it doesn't already exist
    persistent fem;
    if isempty(fem)
        fem = FiniteElementsModel(params);
    end
    FEM = fem;
end
