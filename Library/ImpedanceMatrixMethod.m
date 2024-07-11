classdef ImpedanceMatrixMethod
    properties (Access = public)
        steady = struct('delta_t', [], 'Ut', dictionary());
    end

    properties (Access = private)
        % DOF of harmonic (uh) nodes after removing the ground nodes
        uh = [];
        % All excitation frequencies (w)
        w = [];
    end

    methods (Access = public)
        function [obj] = ImpedanceMatrixMethod(Modal)
            disp('---- Impedance matrix method');
            obj.w = Modal.w;
            obj.uh = Modal.uh;
        end

        function [obj] = computeForcedHarmonicResponse(obj, params, FEM)
            obj = obj.impedanceSteadyStateForcedHarmonicResponse(params, FEM);
        end
    end

    methods (Access = private)
        function [obj] = impedanceSteadyStateForcedHarmonicResponse(obj, params, FEM)
            % Compute the displacement response of harmonic excitation
            % using impedance matrix method
            disp('-- Computing steady-state response using impedance matrix...');

            F = obj.applyExcitationForce(FEM, params.nodes_harmonic, params.Fc, params.Fs);
            obj.steady.delta_t = obj.determineDurationTime(max(obj.w), 0, 2*pi/obj.w(2));
            
            disp('---- Finding response for specified nodes');
            for i = 1:length(params.nodes_harmonic)
                Ut = struct('x', zeros(length(obj.steady.delta_t), length(obj.w)), ...
                            'A', zeros(length(obj.w), 1), 'q', zeros(length(obj.w), 1));
                for j = 1:length(obj.w)
                    Hwj = obj.computeImpedanceFrequencyFunction(FEM, obj.w(j));
                    xt = reshape(Hwj \ F, int32(length(F)/2.0), []);
                    Aj = xt(obj.uh(i), 1);
                    Bj = xt(obj.uh(i), 2);
                    [Ut.x(:,j), Ut.A(j), Ut.q(j)] = obj.computeImpedanceDisplacement(obj.w(j), Aj, Bj);
                end
                obj.steady.Ut(sprintf("%d-%d", params.nodes_harmonic(i), obj.uh(i))) = Ut;
            end
        end

        function [F] = applyExcitationForce(obj, FEM, nodes, Fc, Fs)
            % Apply the excitation force to the specified degrees of freedom
            disp('---- Applying the excitation force');
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

        function [x, A, q] = computeImpedanceDisplacement(obj, wj, Aj, Bj)
            % Compute the displacement response (x), amplitude (A), and phase (q) for excitation wj
            x = Aj * cos(wj * obj.steady.delta_t) + Bj * sin(wj * obj.steady.delta_t);
            A = sqrt(Aj^2 + Bj^2);
            q = abs(atan2(Bj, Aj));
        end

        function [Hwj] = computeImpedanceFrequencyFunction(obj, FEM, wj)
            % Compute frequency response function for the excitation frequency wj
            % Block matrices of the impedance matrix Hwj = [H11, H12; H21, H22]
            H11 =  FEM.Kff - FEM.Mff * wj^2;
            H12 =  FEM.Cff * wj;
            H21 = -FEM.Cff * wj;
            H22 =  FEM.Kff - FEM.Mff * wj^2;
            Hwj = [H11, H12; H21, H22];
        end
    end
end
