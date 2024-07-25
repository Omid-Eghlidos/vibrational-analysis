function [CB] = craigBamptonMethod(FEM)
    fprintf('%s Craig-Bampton method\n', repmat('-', 1, 4));
    CB = getCraigBamptonInstance(FEM);
    % Selecting and adding the substructure
    CB = CB.addSubstructures(FEM);
    % Perform the Craig-Bampton reduction
    CB = CB.performCraigBamptonReduction();
end


function [Craig_Bampton] = getCraigBamptonInstance(FEM)
    % Create an instance of CraigBampton class only if it doesn't already exist
    fprintf('%s Creating an instance of the Craig Bampton class\n', repmat('-', 1, 6));
    persistent craig_bampton;
    if isempty(craig_bampton)
        craig_bampton = CraigBampton(FEM);
    else
        fprintf('%s Already created - skip\n', repmat('-', 1, 8))
    end
    Craig_Bampton = craig_bampton;
end
