classdef ImpedanceMatrix
    properties (Access = public)
        steady = struct('delta_t', [], 'Ut', dictionary());
        % All excitation frequencies (w)
        w = [];
    end

    properties (Access = private)
        % DOF of harmonic (uh) nodes after removing the ground nodes
        uh = [];
    end

    methods (Access = public)
        function [obj] = ImpedanceMatrix(Modal)
            obj.w = Modal.w;
            obj.uh = Modal.uh;
        end

        function [obj] = computeForcedHarmonicResponse(obj, params, FEM)
            % Initializing a parallel pool (used for computing impedance matrix response)
            if isempty(gcp('nocreate'))
                parpool;
            end
            obj = obj.impedanceSteadyStateForcedHarmonicResponse(params, FEM);
        end
    end

    methods (Access = private)
        function [obj] = impedanceSteadyStateForcedHarmonicResponse(obj, params, FEM)
            % Compute the displacement response of harmonic excitation
            % using impedance matrix method
            fprintf('%s Computing steady-state damped response using impedance matrix\n', repmat('-', 1, 6));
            obj.steady.delta_t = obj.determineDurationTime(max(obj.w), 0, 2*pi/obj.w(2));
            F = obj.applyExcitationForce(FEM, params.nodes_harmonic, params.Fc, params.Fs);
            
            for i = 1:length(params.nodes_harmonic)
                Ut = struct('x', [], 'A', [], 'q', []);
                [Ut.x, Ut.A, Ut.q] = obj.computeImpedanceDisplacement(FEM, F, obj.uh(i));      
                obj.steady.Ut(sprintf("%d-%d", params.nodes_harmonic(i), obj.uh(i))) = Ut;
            end
        end

        function [F] = applyExcitationForce(obj, FEM, nodes, Fc, Fs)
            % Apply the excitation force to the specified degrees of freedom
            % Cosine forces
            Fcos = zeros(length(FEM.Mff), 1);
            if ~isempty(Fc)
                [~, nnFc] = ismember(Fc(1,1), nodes(:,1));
                Fcos(obj.uh(nnFc), 1) = Fc(1,3);
            end
            % Sine forces
            Fsin = zeros(length(FEM.Mff), 1);
            if ~isempty(Fs)
                [~, nnFs] = ismember(Fs(1,1), nodes(:,1));
                Fsin(obj.uh(nnFs), 1) = Fs(1,3);
            end
            % Total forces
            F = [Fcos; Fsin];
        end

        function [delta_t] = determineDurationTime(obj, w, t_start, t_end)
            % Determine the duration of the excitation force
            % Number of points over each delta t to interpolate
            ndt = 72;
            T = 2*pi / w;
            dt = T / ndt;
            delta_t = (t_start:dt:t_end)';
        end

        function [x, A, q] = computeImpedanceDisplacement(obj, FEM, F, dof)
            % Compute the displacement response (x), amplitude (A), and phase (q) for excitation wj
            x = zeros(length(obj.steady.delta_t), length(obj.w));
            A = zeros(length(obj.w), 1);
            q = zeros(length(obj.w), 1);
            % Start parallel worker to compute response for each wj
            parfor j = 1:length(obj.w)
                % Compute frequency response function for the excitation frequency wj
                % Block matrices of the impedance matrix Hwj = [H11, H12; H21, H22]
                H11 =  FEM.Kff - FEM.Mff * obj.w(j)^2;
                H12 =  FEM.Cff * obj.w(j);
                H21 = -FEM.Cff * obj.w(j);
                H22 =  FEM.Kff - FEM.Mff * obj.w(j)^2;
                Hwj = [H11, H12; H21, H22];
                xt = reshape(Hwj \ F, int32(length(F)/2.0), []);
                Aj = xt(dof, 1);
                Bj = xt(dof, 2);
                x(:,j) = Aj * cos(obj.w(j) * obj.steady.delta_t) + Bj * sin(obj.w(j) * obj.steady.delta_t);
                A(j) = sqrt(Aj^2 + Bj^2);
                q(j) = abs(atan2(Bj, Aj));
            end
        end
    end
end
