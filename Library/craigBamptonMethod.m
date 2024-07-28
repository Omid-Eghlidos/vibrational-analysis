function [CB] = craigBamptonMethod(params, FEM)
    fprintf('%s Craig-Bampton method\n', repmat('-', 1, 4));
    CB = getCraigBamptonInstance(params);
    % Selecting and adding the substructure
    CB = CB.addSubstructures(FEM);
    % Perform the Craig-Bampton reduction
    CB = CB.performCraigBamptonReduction();
end


function [Craig_Bampton] = getCraigBamptonInstance(params)
    % Create an instance of CraigBampton class only if it doesn't already exist
    persistent craig_bampton;
    if isempty(craig_bampton)
        craig_bampton = CraigBampton(params);
    end
    Craig_Bampton = craig_bampton;
end
