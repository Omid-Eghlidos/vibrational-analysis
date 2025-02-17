classdef ModalExpansion
    properties (Access = public)
        % Eigenvalues and eigenvectores (not sorted)
        Phi = []; Lambda = [];
        % Sorted eigenvalues (natural frequencies) and eigenvectors
        phi = []; fn = [];
        % fn in the given frequency range for damped forced responses
        fn_harmonic = [];
        % DOF of harmonic (uh) nodes after removing the ground nodes
        uh = []; 
        % All the excitation frequencies (w)
        w = [];
        % Free and ground nodes low frequency vibration shapes
        xg = []; yg = []; zg = [];
        xf = []; yf = []; zf = [];
        % Damped steady-state, transient, and free responses
        steady    = struct('delta_t', [], 'Ut', dictionary(), 'Vt', dictionary(), 'At', dictionary(), 'Rt', dictionary());
        transient = struct('delta_t', [], 'Ut', dictionary(), 'Vt', dictionary());
        forced    = struct('delta_t', [], 'Ut', dictionary(), 'Vt', dictionary());
        free      = struct('delta_t', [], 'Ut', dictionary(), 'Vt', dictionary());
    end

    properties (Access = private)
        % wn in the given frequency range for forced and free responses
        wn_harmonic = []; wn_free = [];
        % Significant modes based on the participation factor
        modes = [];
        % DoF of transient (ut) nodes after removing the ground nodes
        ut = [];
        % Modal mass, damping, and stiffness matrices
        MM = []; KK = []; CC = []; FF = [];
        % Damping ratios for natural modes and excitations
        zeta = struct('modes', [], 'excitations', []);
    end

    methods (Access = public)
        function [obj] = ModalExpansion(params, FEM)
            obj = obj.findAndSeparateModes(params, FEM);
            obj = obj.separateExcitedNodes(params, FEM);
            obj = obj.findSignificantModes(params, FEM);
            obj = obj.lowFrequencyVibrationShapes(params, FEM);
        end

        function [obj] = computeForcedHarmonicResponse(obj, params, FEM)
            obj = obj.determineExcitationFrequenciesAndDampingRatios(params);
            obj = obj.dampedSteadyStateForcedHarmonicResponse(params, FEM);
            obj = obj.dampedTransientForcedHarmonicResponse(params, FEM);
            obj = obj.dampedForcedHarmonicResponse(params, FEM);
        end

        function [obj] = computeDampedFreeResponse(obj, params)
            obj = obj.dampedFreeResponse(params);
        end
    end

    methods (Access = private)
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
            fprintf('%s Computing steady-state damped response\n', repmat('-', 1, 6));
            obj.steady.delta_t = obj.determineExcitationTime(max(obj.w), 0, 2*pi/obj.w(2));
            F = obj.applyExcitationForce(params.nodes_harmonic, params.Fc);
            modal = obj.computeModalMatrices(FEM, F, length(obj.wn_harmonic));

            for i = 1:length(params.nodes_harmonic)
                [Ut, Vt, At, Rt] = obj.initializeResponseMatrices(length(obj.steady.delta_t), length(obj.w));
                for j = 1:length(obj.w)
                    [xt, Aj, Bj] = obj.computeResponseCoefficients(modal, obj.w(j), obj.uh(i));
                    [Ut.x(:,j), Ut.A(j), Ut.q(j)] = obj.computeDampedSteadyStateDisplacement(obj.w(j), Aj, Bj, obj.steady.delta_t);
                    [Vt.x(:,j), Vt.A(j), Vt.q(j)] = obj.computeDampedSteadyStateVelocity(obj.w(j), Aj, Bj, obj.steady.delta_t);
                    [At.x(:,j), At.A(j), At.q(j)] = obj.computeDampedSteadyStateAcceleration(obj.w(j), Aj, Bj, obj.steady.delta_t);
                    [Rt.x(:,j), Rt.A(j), Rt.q(j)] = obj.computeReactionForce(FEM, obj.w(j), xt, obj.steady.delta_t);
                end
                obj.steady.Ut(sprintf("%d-%d", params.nodes_harmonic(i), obj.uh(i))) = Ut;
                obj.steady.Vt(sprintf("%d-%d", params.nodes_harmonic(i), obj.uh(i))) = Vt;
                obj.steady.At(sprintf("%d-%d", params.nodes_harmonic(i), obj.uh(i))) = At;
                obj.steady.Rt(sprintf("%d-%d", params.nodes_harmonic(i), FEM.ug(3))) = Rt;
            end
        end

        function [obj] = dampedTransientForcedHarmonicResponse(obj, params, FEM)
            % Determine the damped transient response vs. time of a set of
            % points of the structure under harmonic excitations for sepcified
            % frequencies and zero initial conditions using only mode 1
            fprintf('%s Computing transient damped response\n', repmat('-', 1, 6));
            obj.transient.delta_t = obj.determineExcitationTime(max(obj.w), 0, params.time0);
            F = obj.applyExcitationForce(params.nodes_harmonic, params.Fc);
            modal = obj.computeModalMatrices(FEM, F, length(obj.wn_harmonic));

            for i = 1:length(params.nodes_transient)
                [Ut, Vt] = obj.initializeResponseMatrices(length(obj.transient.delta_t), length(obj.w));
                for j = 1:length(obj.w)
                    zwj = obj.zeta.excitations(j) * obj.w(j);
                    wdj = (1 - 5*obj.zeta.excitations(j)) * obj.w(j);
                    [~, Aj, Bj] = obj.computeResponseCoefficients(modal, zwj, obj.uh(i));
                    [Ut.x(:,j), Ut.A(j)] = obj.computeDampedTransientDisplacement(zwj, wdj, Aj, Bj, obj.transient.delta_t);
                    [Vt.x(:,j), Vt.A(j)] = obj.computeDampedTransientVelocity(Ut.x(:,j), zwj, wdj, Aj, Bj, obj.transient.delta_t);
                end
                obj.transient.Ut(sprintf("%d-%d", params.nodes_transient(i), obj.ut(i))) = Ut;
                obj.transient.Vt(sprintf("%d-%d", params.nodes_transient(i), obj.ut(i))) = Vt;
            end
        end

        function [obj] = dampedForcedHarmonicResponse(obj, params, FEM)
            % Compute damped force response (transient + steady-state) vs. time
            % of a set of points of the structure under harmonic excitations
            % zero initial conditions
            fprintf('%s Computing total forced damped response\n', repmat('-', 1, 6));
            obj.forced.delta_t = obj.determineExcitationTime(max(obj.w), 0, params.time0);
            F = obj.applyExcitationForce(params.nodes_harmonic, params.Fc);
            modal = obj.computeModalMatrices(FEM, F, length(obj.wn_harmonic));

            for i = 1:length(params.nodes_transient)
                [Ut, Vt] = obj.initializeResponseMatrices(length(obj.forced.delta_t), length(obj.w));
                for j = 1:length(obj.w)
                    zwj = obj.zeta.excitations(j) * obj.w(j);
                    wdj = (1 - 5*obj.zeta.excitations(j)) * obj.w(j);
                    % Coefficients for steady-state and transient responses
                    Aj = struct(); Bj = struct();
                    [~, Aj.ss, Bj.ss] = obj.computeResponseCoefficients(modal, obj.w(j), obj.uh(i));
                    [~, Aj.tr, Bj.tr] = obj.computeResponseCoefficients(modal, zwj, obj.uh(i));
                    % Forced harmonic response
                    [Ut.x(:,j), Ut.A(j)] = obj.computeDampedForcedDisplacement(obj.w(j), zwj, wdj, Aj, Bj, obj.forced.delta_t);
                    [Vt.x(:,j), Vt.A(j)] = obj.computeDampedForcedVelocity(Ut.x(:,j), obj.w(j), zwj, wdj, Aj, Bj, obj.transient.delta_t);
                end
                obj.forced.Ut(sprintf("%d-%d", params.nodes_transient(i), obj.ut(i))) = Ut;
                obj.forced.Vt(sprintf("%d-%d", params.nodes_transient(i), obj.ut(i))) = Vt;
            end
        end

        function obj = dampedFreeResponse(obj, params)
            % Determine the damped free response vs. time using the forced
            % response as the initial conditions for displacements and velocities
            fprintf('%s Computing free damped response\n', repmat('-', 1, 6));
            obj.free.delta_t = obj.determineExcitationTime(obj.wn_free(1), params.time0, params.time1);

            for i = 1:length(params.nodes_transient)
                [Ut, Vt] = obj.initializeResponseMatrices(length(obj.free.delta_t), length(obj.wn_free));
                for j = 1:length(obj.w)
                    zwj = obj.zeta.excitations(j) * obj.w(j);
                    wdj = (1 - 5*obj.zeta.excitations(j)) * obj.w(j);
                    % Initial boundary conditions
                    x0 = obj.forced.Ut(sprintf("%d-%d", params.nodes_harmonic(i), obj.ut(i))).x(end,j);
                    v0 = obj.forced.Vt(sprintf("%d-%d", params.nodes_harmonic(i), obj.ut(i))).x(end,j);
                    Aj = x0;
                    Bj = (v0 + zwj * x0) / wdj;
                    [Ut.x(:,j), Ut.A(j)] = obj.computeDampedFreeDisplacement(zwj, wdj, Aj, Bj, obj.free.delta_t);
                    [Vt.x(:,j), Vt.A(j)] = obj.computeDampedFreeVelocity(Ut.x(:,j), zwj, wdj, Aj, Bj, obj.free.delta_t);
                end
                obj.free.Ut(sprintf("%d-%d", params.nodes_transient(i), obj.ut(i))) = Ut;
                obj.free.Vt(sprintf("%d-%d", params.nodes_transient(i), obj.ut(i))) = Vt;
            end
        end

        %% Helper methods
        function [obj] = findAndSeparateModes(obj, params, FEM)
            % Find the vibrational modes and divide them into free, harmonic,
            % and transient modes in the given ranges
            fprintf('%s Finding modes and mode shapes\n', repmat('-', 1, 6));
            % Finding and sort the eigenvalues and eigenvectors for free nodes
            [obj.Phi, obj.Lambda] = eig(FEM.Kff, FEM.Mff);
            [w2, idx] = sort(diag(obj.Lambda));
            obj.phi = obj.Phi(:, idx);

            % Find all the natural frequencies (fn)
            obj.fn = sqrt(w2) / (2*pi);
            % Find wn and fn in the given free range
            fn_free = obj.fn(obj.fn >= min(params.free_range) & obj.fn <= max(params.free_range));
            obj.wn_free = 2*pi * fn_free;
            % Find wn and fn in the given harmonic range
            obj.fn_harmonic = obj.fn(obj.fn >= params.range_harmonic(1) & obj.fn <= params.range_harmonic(2));
            obj.wn_harmonic = 2*pi * obj.fn_harmonic;
        end

        function [obj] = separateExcitedNodes(obj, params, FEM)
            % Find the new DoFs after removal of the ground nodes
            % Harmonic nodes
            [~, nnh] = ismember(params.nodes_harmonic(:,1), FEM.nnf);
            obj.uh = (nnh - 1) * 6 + params.nodes_harmonic(:,2);
            % Transient nodes
            [~, nnt] = ismember(params.nodes_transient(:,1), FEM.nnf);
            obj.ut = (nnt - 1) * 6 + params.nodes_transient(:,2);
        end

        function [obj] = findSignificantModes(obj, params, FEM)
            % Find the significant modes using participation factor of each mode
            fprintf('%s Finding significant modes\n', repmat('-', 1, 6));
            [~, Meff] = obj.computeParticipationFactorAndEffectiveMass(params, FEM);
            % Sort effective masses in descending order and calculate cumulative sum
            [Meff_sorted, sorted_indices] = sort(Meff, 'descend');
            Meff_cumulative = cumsum(Meff_sorted);
            % Total mass of the system
            Mtot = sum(diag(FEM.Mff));
            % Determine the threshold percentage of the total mass (e.g., 90% of the total mass)
            tol = 0.9;
            Mthreshold = tol * Mtot;
            % Find the number of modes required to reach the threshold
            modes_to_include = sorted_indices(1:find(Meff_cumulative >= Mthreshold, 1));
            % If no significant modes meet the threshold, use a fallback number of modes (e.g., 10)
            if isempty(modes_to_include)                
                modes_to_include = length(obj.wn_harmonic);
                fprintf('%s No modes meet the threshold, fall back to the first %d modes\n', repmat('-', 1, 8), modes_to_include);
                obj.modes = (1:modes_to_include);
            else
                obj.modes = sorted_indices(1:modes_to_include);
            end
            
        end

        function [gamma, Meff] = computeParticipationFactorAndEffectiveMass(obj, params, FEM)
            obj = obj.checkModeVectorNormalization(FEM);
            % Create the excitation direction factor (r)
            r = zeros(length(FEM.Mff), 1);
            r(obj.uh(1)) = params.Fc(3);
            % Compute participation factors and effective masses
            gamma = zeros(length(FEM.Mff), 1);
            Meff = zeros(length(FEM.Mff), 1);
            for j = 1:length(FEM.Mff)
                gamma(j) = (obj.phi(:, j)' * FEM.Mff * r) / (obj.phi(:, j)' * FEM.Mff * obj.phi(:, j));
                Meff(j) = (obj.phi(:, j)' * FEM.Mff * r)^2 / (obj.phi(:, j)' * FEM.Mff * obj.phi(:, j));
            end
        end

        function [obj] = checkModeVectorNormalization(obj, FEM)
            % Normalize the mode vectors with respect to mass matrix if needed
            normalized = true;
            for j = 1:length(FEM.Mff)
                normalization_value = round(obj.phi(:,j)' * FEM.Mff * obj.phi(:,j), 6);
                if abs(normalization_value - 1) > 1e-6
                    normalized = false;
                    break;
                end
            end

            if ~normalized
                obj = obj.normalizeModeVectors(FEM);
            end
        end

        function [obj] = normalizeModeVectors(obj, FEM)
            % Normalize mode vectors with respect to the mass matrix
            for j = 1:length(FEM.Mff)
                obj.phi(:,j) = obj.phi(:,j) / sqrt(obj.phi(:,j)' * FEM.Mff * obj.phi(:,j));
            end
        end

        function [obj] = determineExcitationFrequenciesAndDampingRatios(obj, params)
            % Determine excitation frequencies by interpolating the given range
            % To reduce computational cost, interpolate regions around
            % each natural frequency (wn) with more points than regions further away
            % Number of points to interpolate the half-power bandwidth (even number)
            nbw = 6;
            % Fine resolution frequencies (wf) for all the bands
            wf = zeros(nbw * length(obj.wn_harmonic), 1);
            obj.zeta.modes = zeros(length(obj.wn_harmonic), 1);
            for mode = 1:length(obj.wn_harmonic)
                % Compute damping ratio for each mode using Rayleigh
                obj.zeta.modes(mode) = params.alpha / (2*obj.wn_harmonic(mode)) + params.beta * (obj.wn_harmonic(mode)/2);
                % Half-bandwidth points on both sides of the peak (wn)
                w1 = obj.wn_harmonic(mode) * (1 - obj.zeta.modes(mode));
                w2 = obj.wn_harmonic(mode) * (1 + obj.zeta.modes(mode));
                wf(nbw*(mode-1)+1:nbw*(mode)+1, 1) = linspace(w1, w2, nbw+1);
            end

            % Determine number of frequencies for coarse resolution
            num_wc = params.num_w - nbw * length(obj.wn_harmonic);
            if num_wc < 0
                num_wc = 2*length(wf);
            end
            % Interpolate the rest with less points - Coarse resolution frequencies (wc)
            wc = linspace(0, 2*pi*params.range_harmonic(end), num_wc)';
            % Total excitation frequencies (w)
            obj.w = unique(sort([wc; wf; obj.wn_harmonic]), 'stable');
            % Interpolate zeta for the excitation frquencies
            obj.zeta.excitations = interp1(2*pi*params.zetas(:,1), params.zetas(:,2), obj.w, 'linear', 'extrap');
        end

        function [F] = applyExcitationForce(obj, nodes, Force)
            % Apply the excitation force to the specified degrees of freedom
            F = zeros(length(obj.phi), 1);
            [~, nnF0] = ismember(Force(1,1), nodes(:,1));
            F(obj.uh(nnF0), 1) = Force(1,3);
        end

        function [delta_t] = determineExcitationTime(obj, w, t_start, t_end)
            % Determine the duration of the excitation force
            % Number of points to interpolate time (multiply of 36)
            ndt = 72;
            T = 2 * pi / w;
            dt = T / ndt;
            delta_t = (t_start:dt:t_end)';
        end

        function [modal] = computeModalMatrices(obj, FEM, F, modes)
            % Compute modal mass, damp, stiffness matrices for the given DoFs
            modal = struct();
            modal.M = obj.phi(:,1:modes)' * FEM.Mff * obj.phi(:,1:modes);
            modal.C = obj.phi(:,1:modes)' * FEM.Cff * obj.phi(:,1:modes);
            modal.K = obj.phi(:,1:modes)' * FEM.Kff * obj.phi(:,1:modes);
            modal.F = obj.phi(:,1:modes)' * F;
        end

        function [Ut, Vt, varargout] = initializeResponseMatrices(obj, m, n)
            % Initialize matrices for displacement (Ut), velocity (Vt),
            % acceleration (At), and reaction force (Rt) for each node
            Ut = struct('x', zeros(m, n), 'A', zeros(n, 1));
            Vt = struct('x', zeros(m, n), 'A', zeros(n, 1));
            if nargout > 2
                Ut.q = zeros(n, 1);
                Vt.q = zeros(n, 1);
                At = struct('x', zeros(m, n), 'A', zeros(n, 1), 'q', zeros(n, 1));
                Rt = struct('x', zeros(m, n), 'A', zeros(n, 1), 'q', zeros(n, 1));
                varargout{1} = At;
                varargout{2} = Rt;
            end
        end

        function [xt, Aj, Bj] = computeResponseCoefficients(obj, modal, wj, dof)
            % Compute frequency response function for excitation frequency wj
            Hwj = -wj^2 * modal.M + 1i * wj * modal.C + modal.K;
            xt = obj.phi(:,1:length(obj.wn_harmonic)) * (Hwj \ modal.F);
            Aj = real(xt(dof));
            Bj = imag(xt(dof));
        end

        function [x, A, q] = computeDampedSteadyStateDisplacement(obj, wj, Aj, Bj, t)
            % Compute steady-state displacement response (x), amplitude (A), 
            % and phase (q) for excitation wj
            x = Aj * cos(wj * t) + Bj * sin(wj * t);
            A = sqrt(Aj^2 + Bj^2);
            q = abs(atan2(Bj, Aj));
        end

        function [x, A, q] = computeDampedSteadyStateVelocity(obj, wj, Aj, Bj, t)
            % Compute steady-state velocity response (x), amplitude (A), 
            % and phase (q) for excitation wj
            x = -wj * Aj * sin(wj * t) + wj * Bj * cos(wj * t);
            A = wj * sqrt(Aj^2 + Bj^2);
            q = abs(atan2(wj * Bj, wj * Aj));
        end

        function [x, A, q] = computeDampedSteadyStateAcceleration(obj, wj, Aj, Bj, t)
            % Compute steady-state acceleration response (x), amplitude (A), 
            % and phase (q) for excitation wj
            x = -wj^2 * Aj * cos(wj * t) - wj^2 * Bj * sin(wj * t);
            A = wj^2 * sqrt(Aj^2 + Bj^2);
            q = abs(atan2(wj^2 * Bj, wj^2 * Aj));
        end

        function [x, A, q] = computeReactionForce(obj, FEM, wj, xt, t)
            % Compute the reaction forces' response (x), amplitude (A), 
            % and phase (q) for excitation wj
            dt = t(2) - t(1);
            vt = gradient(xt, dt);
            at = gradient(vt, dt);
            r = FEM.Mgf * at + FEM. Cgf * vt + FEM.Kgf * xt;
            % Get the force on the first ground node along its 3rd DoF
            Fc = real(r(FEM.ug(3)));
            Fs = imag(r(FEM.ug(3)));
            x = Fc * cos(wj * t) + Fs * sin(wj * t);
            A = sqrt(Fc^2 + Fs^2);
            q = abs(atan2(Fs, Fc));
        end

        function [x, A] = computeDampedTransientDisplacement(obj, zwj, wdj, Aj, Bj, t)
            % Compute transient displacement response (x) and amplitude (A) for excitation wj
            x = exp(-zwj * t) .* (Aj * cos(wdj * t) + Bj * sin(wdj * t));
            A = sqrt(Aj^2 + Bj^2);
        end

        function [x, A] = computeDampedTransientVelocity(obj, Ut, zwj, wdj, Aj, Bj, t)
            % Compute transient velocity response (x), amplitude (A), and phase (q) for excitation wj
            x = -zwj * Ut + wdj * exp(-zwj * t) .* (-Aj * sin(wdj * t) + Bj * cos(wdj * t));
            A = wdj * sqrt(Aj^2 + Bj^2);
        end

        function [x, A] = computeDampedForcedDisplacement(obj, wj, zwj, wdj, Aj, Bj, t)
            % Compute forced (steady-state + transient) displacement response (x)
            % and amplitude (A) for excitation wj
            % Steady-state response
            [xss, ~, ~] = obj.computeDampedSteadyStateDisplacement(wj, Aj.ss, Bj.ss, t);
            A = sqrt(Aj.ss^2 + Bj.ss^2);
            % Transient response
            [xtr, ~] = obj.computeDampedTransientDisplacement(zwj, wdj, Aj.tr, Bj.tr, t);
            A = A + sqrt(Aj.tr^2 + Bj.tr^2);
            x = xtr + xss;
        end

        function [x, A] = computeDampedForcedVelocity(obj, Ut, wj, zwj, wdj, Aj, Bj, t)
            % Compute forced (steady-state + transient) displacement response (x)
            % and amplitude (A) for excitation wj
            % Steady-state response
            [xss, ~, ~] = obj.computeDampedSteadyStateVelocity(wj, Aj.ss, Bj.ss, t);
            A = wj * sqrt(Aj.ss^2 + Bj.ss^2);
            % Transient response
            [xtr, ~] = obj.computeDampedTransientVelocity(Ut, zwj, wdj, Aj.tr, Bj.tr, t);
            A = A + wj * sqrt(Aj.tr^2 + Bj.tr^2);
            x = xtr + xss;
        end

        function [x, A] = computeDampedFreeDisplacement(obj, zwj, wdj, Aj, Bj, t)
            % Compute free displacement response (x) and amplitude (A) for excitation wj
            x = exp(-zwj * t) .* (Aj * cos(wdj * t) + Bj * sin(wdj * t));
            A = sqrt(Aj^2 + Bj^2);
        end

        function [x, A] = computeDampedFreeVelocity(obj, Ut, zwj, wdj, Aj, Bj, t)
            % Compute free velocity response (x) and amplitude (A) for excitation wj
            x = -zwj * Ut + wdj * exp(-zwj * t) .* (-Aj * sin(wdj * t) + Bj * cos(wdj * t));
            A = wdj * sqrt(Aj^2 + Bj^2);
        end

    end
end
