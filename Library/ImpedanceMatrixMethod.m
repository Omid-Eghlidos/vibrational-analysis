function [Impedance] = impedanceMatrixMethod(params, FEM, Modal)
    fprintf('%s Impedance matrix method\n', repmat('-', 1, 4));
    Impedance = getImpedanceMatrixInstance(Modal);
    % Compute the damped steady-state response using the impedance matrix method
    Impedance = Impedance.computeForcedHarmonicResponse(params, FEM);
end


function [Impedance] = getImpedanceMatrixInstance(Modal)
    % Create an instance of ImpedanceMatrix class only if it doesn't already exist
    persistent impedance;
    if isempty(impedance)
        impedance = ImpedanceMatrix(Modal);
    end
    Impedance = impedance;
end
