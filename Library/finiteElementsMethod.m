function [FEM] = finiteElementsMethod(params)
    fprintf('%s Finite elements method\n', repmat('-', 1, 4));
    FEM = getFiniteElementsModelInstance(params);
end


function [FEM] = getFiniteElementsModelInstance(params)
    % Create an instance of FiniteElementsModel class only if it doesn't already exist
    fprintf('%s Creating an instance of the Finite Elements Model class\n', repmat('-', 1, 6));
    persistent fem;
    if isempty(fem)
        fem = FiniteElementsModel(params);
    else
        fprintf('%s Already created - skip\n', repmat('-', 1, 8))
    end
    FEM = fem;
end
