function [Modal] = modalExpansionMethod(params, FEM)
    fprintf('%s Modal expansion method\n', repmat('-', 1, 4));
    Modal = getModalExpansionInstance(params, FEM);
    % Compute the damped steady-state and transient responses for harmonic excitations
    Modal = Modal.computeForcedHarmonicResponse(params, FEM);
    % Compute the damped free response with final value of forced response as its initial condition
    Modal = Modal.computeDampedFreeResponse(params);
end


function [Modal] = getModalExpansionInstance(params, FEM)
    % Create an instance of ModalExpansion class only if it doesn't already exist
    persistent modal;
    if isempty(modal)
        modal = ModalExpansion(params, FEM);
    end
    Modal = modal;
end
