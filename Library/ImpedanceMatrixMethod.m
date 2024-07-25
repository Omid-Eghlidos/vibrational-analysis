function [Impedance] = impedanceMatrixMethod(params, FEM, Modal)
    fprintf('%s Impedance matrix method\n', repmat('-', 1, 4));
    Impedance = getImpedanceMatrixInstance(Modal);
    % Compute the damped steady-state response using the impedance matrix method
    Impedance = Impedance.computeForcedHarmonicResponse(params, FEM);
end


function [Impedance] = getImpedanceMatrixInstance(Modal)
    % Create an instance of ImpedanceMatrix class only if it doesn't already exist
    fprintf('%s Creating an instance of the Impedance Matrix class\n', repmat('-', 1, 6));
    persistent impedance;
    if isempty(impedance)
        impedance = ImpedanceMatrixMethod(Modal);
    else
        fprintf('%s Already created - skip\n', repmat('-', 1, 8))
    end
    Impedance = impedance;
end
