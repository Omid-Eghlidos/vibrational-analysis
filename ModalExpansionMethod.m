classdef ModalExpansionMethod
    properties (Access = public)
        % Sorted eigenvalues (natural frequencies) and eigenvectors
        phi = []; fn = [];
        % fn in the given frequency range for damped forced and free responses
        fn_harmonic = []; fn_free = [];
        % Free and ground nodes low frequency vibration shapes
        xg = []; yg = []; zg = [];
        xf = []; yf = []; zf = [];
        % All the excitation frequencies (w)
        w = [];
        % Damped steady-state, transient, and free responses
        steady    = struct('delta_t', [], 'Ut', dictionary(), 'Vt', dictionary(), 'At', dictionary(), 'Rt', dictionary());
        transient = struct('delta_t', [], 'Ut', dictionary(), 'Vt', dictionary(), 'At', dictionary());
        free      = struct('delta_t', [], 'Ut', dictionary(), 'Vt', dictionary(), 'At', dictionary());
    end

    properties (Access = private)
        % Modal mass, damping, and stiffness matrices
        MM = []; KK = []; CC = []; FF = [];
        % DOF of harmonic (uh) and transient (ut) nodes after removing the ground nodes
        uh = []; ut = [];
        % wn in the given frequency range for forced and free responses
        wn_harmonic = []; wn_free = [];
        % Significant modes based on the participation factor
        significant_modes = [];
    end

    methods (Access = public)
        function [obj] = ModalExpansionMethod(params, FEM)
            disp('---- Modal expansion method');
            obj = obj.findAndSeparateModes(params, FEM);
            obj = obj.separateExcitedNodes(params, FEM);
            obj = obj.findSignificantModes(params, FEM);
            obj = obj.lowFrequencyVibrationShapes(params, FEM);
        end

        function [obj] = computeForcedHarmonicResponse(obj, params, FEM)
            obj = obj.determineExcitationFrequencies(params);
            obj = obj.dampedSteadyStateForcedHarmonicResponse(params, FEM);
            obj = obj.dampedTransientForcedHarmonicResponse(params, FEM);
        end

        function [obj] = computeDampedFreeResponse(obj, params)
            obj = obj.dampedFreeResponse(params);
        end
    end

    methods (Access = private)
        function [obj] = findAndSeparateModes(obj, params, FEM)
            % Find the vibrational modes and divide them into free, harmonic, 
            % and transient modes in the given ranges
            disp('------ Finding modes and mode shapes');
            % Finding and sort the eigenvalues and eigenvectors for free nodes
            [V, D] = eig(FEM.Kff, FEM.Mff);
            [d, idx] = sort(diag(D));
            obj.phi = V(:, idx);

            % Find all the natural frequencies (fn)
            obj.fn = sqrt(d) / (2*pi);
            % Find wn and fn in the given free range
            obj.fn_free = obj.fn(obj.fn >= min(params.free_range) & obj.fn <= max(params.free_range));
            obj.wn_free = 2*pi * obj.fn_free;
            % Find wn and fn in the given harmonic range
            obj.fn_harmonic = obj.fn(obj.fn >= params.range_harm(1) & obj.fn <= params.range_harm(2));
            obj.wn_harmonic = 2*pi * obj.fn_harmonic;
        end

        function [obj] = lowFrequencyVibrationShapes(obj, params, FEM)
            % Find the mode shapes for vibrations of the low frequency modes
            % Coordinates of the grounded nodes
            obj.xg = FEM.X(FEM.nng);
            obj.yg = FEM.Y(FEM.nng);
            obj.zg = zeros(length(FEM.nng), 1);
            % Coordinates of the free nodes
            obj.xf = FEM.X(FEM.nnf);
            obj.yf = FEM.Y(FEM.nnf);
            obj.zf = zeros(int32(length(obj.phi)/6), params.low_modes);

            for mode = 1:params.low_modes
                Tx = obj.phi(1:6:length(obj.phi), mode);
                Ty = obj.phi(2:6:length(obj.phi), mode);
                Tz = obj.phi(3:6:length(obj.phi), mode);
                % Choose the translational dof with the highest maximum value
                if max(abs(Tx)) >= max(max(abs(Ty)), max(abs(Tz)))
                    obj.zf(:,mode) = Tx;
                elseif max(abs(Ty)) >= max(max(abs(Tx)), max(abs(Tz)))
                    obj.zf(:,mode) = Ty;
                elseif max(abs(Tz)) >= max(max(abs(Tx)), max(abs(Ty)))
                    obj.zf(:,mode) = Tz;
                end
            end
        end

        function [obj] = dampedSteadyStateForcedHarmonicResponse(obj, params, FEM)
            % Computes damped steady state response and its corresponding
            % amplitudes and phases, as well as the reaction force at
            % the grounded nodes due to harmonic excitations
            disp('-- Computing damped steady-state response...');

            obj.steady.delta_t = obj.determineDurationTime(max(obj.w), 0, params.time0);
            F = obj.applyExcitationForce(params.nodes_harm, params.Fc0);
            obj = computeModalMatrices(obj, FEM, F, length(obj.wn_harmonic));

            disp('---- Finding response for specified nodes and frequencies');
            for i = 1:length(params.nodes_harm)
                [Ut, Vt, At, Rt] = obj.initializeResponseMatrices(length(obj.steady.delta_t), length(obj.w));
                for j = 1:length(obj.w)
                    Hwj = obj.computeFrequencyFunction(obj.w(j));
                    xt = obj.phi(:,1:length(obj.wn_harmonic)) * Hwj * obj.FF;
                    Aj = real(xt(obj.uh(i)));
                    Bj = imag(xt(obj.uh(i)));
                    [Ut.x(:,j), Ut.A(j), Ut.q(j)] = obj.computeDampedSteadyStateDisplacement(obj.w(j), Aj, Bj);
                    [Vt.x(:,j), Vt.A(j), Vt.q(j)] = obj.computeDampedSteadyStateVelocity(obj.w(j), Aj, Bj);
                    [At.x(:,j), At.A(j), At.q(j)] = obj.computeDampedSteadyStateAcceleration(obj.w(j), Aj, Bj);
                    [Rt.x(:,j), Rt.A(j), Rt.q(j)] = obj.computeReactionForce(FEM, obj.w(j), xt);
                end
                obj.steady.Ut(sprintf("%d-%d", params.nodes_harm(i), obj.uh(i))) = Ut;
                obj.steady.Vt(sprintf("%d-%d", params.nodes_harm(i), obj.uh(i))) = Vt;
                obj.steady.At(sprintf("%d-%d", params.nodes_harm(i), obj.uh(i))) = At;
                obj.steady.Rt(sprintf("%d-%d", params.nodes_harm(i), FEM.ug(3))) = Rt;
            end
        end

        function [obj] = dampedTransientForcedHarmonicResponse(obj, params, FEM)
            % Determine the damped transient response vs. time of a set of
            % points of the structure under harmonic excitations for sepcified
            % frequencies and zero initial conditions using only mode 1
            disp('-- Computing damped transient response...');

            zeta = obj.computeDampingRatios(params, "transient");
            obj.transient.delta_t = obj.determineDurationTime(max(obj.w), 0, params.time0);
            F = obj.applyExcitationForce(params.nodes_harm, params.Fc1);
            obj = computeModalMatrices(obj, FEM, F, length(obj.wn_harmonic));

            disp('---- Finding response for specified nodes and frequencies');
            for i = 1:length(params.nodes_trans)
                [Ut, Vt, ~, ~] = obj.initializeResponseMatrices(length(obj.transient.delta_t), length(obj.w));
                for j = 1:length(obj.w)
                    zwj = zeta(j) * obj.w(j);
                    wdj = (1 - 5*zeta(j)) * obj.w(j); 
                    Hwj = obj.computeFrequencyFunction(zwj);
                    xt = obj.phi(:,1:length(obj.wn_harmonic)) * Hwj * obj.FF;
                    Aj = real(xt(obj.uh(i)));
                    Bj = imag(xt(obj.uh(i)));
                    [Ut.x(:,j), Ut.A(j), Ut.q(j)] = obj.computeDampedTransientDisplacement(obj.w(j), zwj, wdj, Aj, Bj);
                    [Vt.x(:,j), Vt.A(j), Vt.q(j)] = obj.computeDampedTransientVelocity(Ut.x(:,j), zwj, wdj, Aj, Bj);
                end
                obj.transient.Ut(sprintf("%d-%d", params.nodes_trans(i), obj.ut(i))) = Ut;
                obj.transient.Vt(sprintf("%d-%d", params.nodes_trans(i), obj.ut(i))) = Vt;
            end
        end

        function obj = dampedFreeResponse(obj, params)
            % Determine the damped free response vs. time using the forced 
            % response as the initial conditions for displacements and velocities 
            disp('-- Computing damped free response...');

            obj.free.delta_t = obj.determineDurationTime(obj.wn_free(1), params.time0, params.time1);

            disp('---- Finding response for specified nodes and frequencies');
            for i = 1:length(params.nodes_trans)
                [Ut, Vt, ~, ~] = obj.initializeResponseMatrices(length(obj.free.delta_t), length(obj.wn_free));
                for j = 1:length(obj.wn_free)
                    zwj = zeta(j) * obj.wn_free(j);
                    wdj = (1 - 5*zeta(j)) * obj.wn_free(j); 
                    % Initial boundary condition
                    x0 = obj.transient.Ut(sprintf("%d-%d", params.nodes_harm(i), obj.ut(i))).x(end,j);
                    v0 = obj.transient.Vt(sprintf("%d-%d", params.nodes_harm(i), obj.ut(i))).x(end,j);
                    Aj = x0;
                    Bj = (v0 + zwj * x0) / wdj;
                    [Ut.x(:,j), Ut.A(j), Ut.q(j)] = obj.computeDampedFreeDisplacement(zwj, wdj, Aj, Bj);
                    [Vt.x(:,j), Vt.A(j), Vt.q(j)] = obj.computeDampedFreeVelocity(Ut.x(:,j), zwj, wdj, Aj, Bj);
                end
                obj.free.Ut(sprintf("%d-%d", params.nodes_trans(i), obj.ut(i))) = Ut;
                obj.free.Vt(sprintf("%d-%d", params.nodes_trans(i), obj.ut(i))) = Vt;
            end
        end

        %% Helper methods
        function [obj] = separateExcitedNodes(obj, params, FEM)
            % Find the new DoFs after removal of the ground nodes
            disp('------ Finding new DoFs of the nodes');
            % Harmonic nodes
            [~, nnh] = ismember(params.nodes_harm(:,1), FEM.nnf);
            obj.uh = (nnh - 1) * 6 + params.nodes_harm(:,2);
            % Transient nodes
            [~, nnt] = ismember(params.nodes_trans(:,1), FEM.nnf);
            obj.ut = (nnt - 1) * 6 + params.nodes_trans(:,2);
        end

        function [obj] = findSignificantModes(obj, params, FEM)
            % Find the significant modes using participation factor of each mode
            disp('------ Finding the significant modes');

            obj.normalizeMassMatrix(FEM);
            % Create the excitation direction factor (r)
            r = zeros(length(FEM.Mff), 1);
            r(obj.uh(1)) = params.Fc0(3);
            % Effective mass vector
            Meff = (obj.phi' * FEM.Mff * obj.phi).^2;
            Meff_cumulative = cumsum(Meff, 2);
            % Total mass of the system
            Mtot = sum(FEM.Mff, 'all');
            % Determine the threshold percentage of the total mass 
            tol = 0.9;
            obj.significant_modes = find(Meff_cumulative(:,end) / Mtot >= tol, 1);
        end

        function [obj] = normalizeMassMatrix(obj, FEM)
            normalized = obj.checkMassMatrixNormalized(FEM);
            if ~normalized
                for i = 1:length(FEM.Mff)
                    obj.phi(:,i) = obj.phi(:,i) / sqrt(obj.phi(:,i)' * FEM.Mff * obj.phi(:,i));
                end
            end
        end

        function [normalized] = checkMassMatrixNormalized(obj, FEM)
            normalized = true;
            % Normalize the mass matrix if needed to obtain the total mass
            for i = 1:length(FEM.Mff)
                normalization_value = round(obj.phi(:,i)' * FEM.Mff * obj.phi(:,i), 6);
                if abs(normalization_value - 1) > 1e-6
                    normalized = false;
                end
            end
        end

        function [obj] = determineExcitationFrequencies(obj, params)
            % Determine excitation frequencies by interpolating the given range
            % To reduce computational cost, interpolate regions around
            % each natural frequency (wn) with more points than regions further away

            zeta = obj.computeDampingRatios(params, "harmonic");
            % Number of points to interpolate in the half-power bandwidth
            nbw = 10;
            % Fine resolution frequencies (wf) for all the bands
            wf = zeros(nbw * length(obj.wn_harmonic), 1);
            for mode = 1:length(obj.wn_harmonic)
                % Half-bandwidth points on both sides of the peak (wn)
                w1 = obj.wn_harmonic(mode) * (1 - zeta(mode));
                w2 = obj.wn_harmonic(mode) * (1 + zeta(mode));
                wf(nbw*(mode-1)+1:nbw*(mode)+1, 1) = linspace(w1, w2, nbw+1);
            end

            % Interpolate the rest with less points - Coarse resolution frequencies (wc)
            wc = linspace(0, 2*pi*params.range_harm(end), length(wf))';
            % Total excitation frequencies (wt) = combine and sort wf, wc, and wn_harmonic
            obj.w = unique(sort([wc; wf; obj.wn_harmonic]), 'stable');
        end

        function [zeta] = computeDampingRatios(obj, params, Type)
            if Type == "harmonic"
                zeta = zeros(length(obj.wn_harmonic), 1);
                for mode = 1:length(obj.wn_harmonic)
                    % Compute damping ratio for each mode using Rayleigh
                    zeta(mode) = params.alpha / (2*obj.wn_harmonic(mode)) + params.beta * (obj.wn_harmonic(mode)/2);
                end
            else
                % Interpolate zeta for frquencies in the given bounds
                zeta = interp1(2*pi*params.damp_ratio_trans(:,1), params.damp_ratio_trans(:,2), obj.w, 'linear', 'extrap');
            end
        end

        function [F] = applyExcitationForce(obj, nodes, Force)
            % Apply the excitation force to the specified degrees of freedom
            disp('---- Applying the excitation force');
            F = zeros(length(obj.phi), 1);
            [~, nnF0] = ismember(Force(1,1), nodes(:,1));
            F(obj.uh(nnF0), 1) = Force(1,3);
        end

        function [delta_t] = determineDurationTime(obj, w, t_start, t_end)
            % Determine the duration of the excitation force
            % Number of points over each delta t to interpolate
            ndt = 72;
            T = 2*pi / w;
            delta_t = (t_start:T/ndt:t_end)';
        end

        function [obj] = computeModalMatrices(obj, FEM, F, modes)
            % Compute modal mass, damp, stiffness matrices for the given DoFs
            disp('---- Computing the modal matrices');
            obj.MM = obj.phi(:,1:modes)' * FEM.Mff * obj.phi(:,1:modes);
            obj.CC = obj.phi(:,1:modes)' * FEM.Cff * obj.phi(:,1:modes);
            obj.KK = obj.phi(:,1:modes)' * FEM.Kff * obj.phi(:,1:modes);
            obj.FF = obj.phi(:,1:modes)' * F;
        end

        function [Hwj] = computeFrequencyFunction(obj, wj)
            % Compute frequency response function for the current excitation frequency
            Hwj = inv(-wj^2 * obj.MM + 1i * wj * obj.CC + obj.KK);
        end

        function [Ut, Vt, At, Rt] = initializeResponseMatrices(obj, m, n)
            % Initialize matrices for displacement (Ut), velocity (Vt), 
            % acceleration (At), and reaction force (Rt) for each node
            Ut = struct('x', zeros(m, n), 'A', zeros(n, 1), 'q', zeros(n, 1));
            Vt = struct('x', zeros(m, n), 'A', zeros(n, 1), 'q', zeros(n, 1));
            At = struct('x', zeros(m, n), 'A', zeros(n, 1), 'q', zeros(n, 1));
            Rt = struct('x', zeros(m, n), 'A', zeros(n, 1), 'q', zeros(n, 1));
        end

        function [x, A, q] = computeDampedSteadyStateDisplacement(obj, wj, Aj, Bj)
            % Compute the displacement response (x), amplitude (A), and phase (q) for excitation wj
            x = Aj * cos(wj * obj.steady.delta_t) + Bj * sin(wj * obj.steady.delta_t);
            A = sqrt(Aj^2 + Bj^2);
            q = abs(atan2(Bj, Aj));
        end

        function [x, A, q] = computeDampedSteadyStateVelocity(obj, wj, Aj, Bj)
            % Compute the velocity response (x), amplitude (A), and phase (q) for excitation wj
            x = -wj * Aj * sin(wj * obj.steady.delta_t) + wj * Bj * cos(wj * obj.steady.delta_t);
            A = wj * sqrt(Aj^2 + Bj^2);
            q = abs(atan2(wj * Bj, wj * Aj));
        end

        function [x, A, q] = computeDampedSteadyStateAcceleration(obj, wj, Aj, Bj)
            % Compute the acceleration response (x), amplitude (A), and phase (q) for excitation wj
            x = -wj^2 * Aj * cos(wj * obj.steady.delta_t) - wj^2 * Bj * sin(wj * obj.steady.delta_t);
            A = wj^2 * sqrt(Aj^2 + Bj^2);
            q = abs(atan2(wj^2 * Bj, wj^2 * Aj));
        end

        function [x, A, q] = computeReactionForce(obj, FEM, wj, xt)
            % Compute the force response (x), amplitude (A), and phase (q) for excitation wj
            dt = obj.steady.delta_t(2) - obj.steady.delta_t(1);
            vt = gradient(xt, dt);
            at = gradient(vt, dt);
            r = FEM.Mgf * at + FEM. Cgf * vt + FEM.Kgf * xt;
            % Get the force on the first ground node along its 3rd DoF
            Fc = real(r(FEM.ug(3)));
            Fs = imag(r(FEM.ug(3)));
            x = Fc * cos(wj * obj.steady.delta_t) + Fs * sin(wj * obj.steady.delta_t);
            A = sqrt(Fc^2 + Fs^2);
            q = abs(atan2(Fs, Fc));
        end

        function [x, A, q] = computeDampedTransientDisplacement(obj, wj, zwj, wdj, Aj, Bj)
            [xss, ~, ~] = obj.computeDampedSteadyStateDisplacement(wj, Aj, Bj);
            % Compute the displacement response (x), amplitude (A), and phase (q) for excitation wj
            x = xss + exp(-zwj * obj.transient.delta_t) .* (Aj * cos(wdj * obj.transient.delta_t) + Bj * sin(wdj * obj.transient.delta_t));
            A = sqrt(Aj^2 + Bj^2);
            q = abs(atan2(Bj, Aj));
        end

        function [x, A, q] = computeDampedTransientVelocity(obj, Ut, zwj, wdj, Aj, Bj)
            % Compute the velocity response (x), amplitude (A), and phase (q) for excitation wj
            x =  -zwj * Ut + wdj * exp(-zwj * obj.transient.delta_t) .* (-Aj * sin(wdj * obj.transient.delta_t) + Bj * cos(wdj * obj.transient.delta_t));
            A = wdj * sqrt(Aj^2 + Bj^2);
            q = abs(atan2(wdj * Bj, wdj * Aj));
        end

        function [x, A, q] = computeDampedFreeDisplacement(obj, zwj, wdj, Aj, Bj)
            % Compute the displacement response (x), amplitude (A), and phase (q) for excitation wj
            x =  exp(-zwj * obj.free.delta_t) .* (Aj * cos(wdj * obj.free.delta_t) + Bj * sin(wdj * obj.free.delta_t));
            A = sqrt(Aj^2 + Bj^2);
            q = abs(atan2(Bj, Aj));
        end

        function [x, A, q] = computeDampedFreeVelocity(obj, Ut, zwj, wdj, Aj, Bj)
            % Compute the velocity response (x), amplitude (A), and phase (q) for excitation wj
            x =  -zwj * Ut + wdj * exp(-zwj * obj.free.delta_t) .* (-Aj * sin(wdj * obj.free.delta_t) + Bj * cos(wdj * obj.free.delta_t));
            A = wdj * sqrt(Aj^2 + Bj^2);
            q = abs(atan2(wdj * Bj, wdj * Aj));
        end

    end
end
