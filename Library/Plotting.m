classdef Plotting
    properties (Access = public)
        folder = dictionary('FEM',  fullfile(pwd, 'Results/FEM'), ...
                            'Modal', fullfile(pwd, 'Results/Modal'), ...
                            'Impedance', fullfile(pwd, 'Results/Impedance'), ...
                            'Craig_Bampton', fullfile(pwd, 'Results/Craig_Bampton'));
        % Plotting modes
        modes2plot = [];
        % All the excitation (f) and harmonic (fn_harmonic) frequencies
        f = []; fn_harmonic = [];
        % Units for response and amplitude of each type of results
        units = dictionary("Displacement", "(mm)", "Velocity", "(mm/s)", ...
                           "Acceleration", "(mm/s$^2$)", "Reaction", "")
    end

    methods (Access = public)
        function [obj] = Plotting()
            vals = values(obj.folder);
            for i = 1:length(keys(obj.folder))
                mkdir(vals(i));
            end
        end

        function finiteElementsModel(obj, params, FEM)
            if params.model
                fprintf('%s Plotting FEM mesh\n', repmat('-', 1, 6));
                obj.finiteElementsMesh(FEM);
                obj.finiteElementsMeshNodes(FEM);
            end
        end

        function modalExpansionMethod(obj, params, Modal)
            fprintf('%s Plotting modal expansion results\n', repmat('-', 1, 6));
            obj.modes2plot = params.modes2plot;
            obj.fn_harmonic = Modal.fn_harmonic;
            obj.f = Modal.w / 2 / pi;
            if params.frequencies
                fprintf('%s Frequency distribution\n', repmat('-', 1, 8));
                obj.frequencyDistribution(Modal, obj.folder('Modal'));
            end
            if params.low_modes
                fprintf('%s Low frequency mode shapes\n', repmat('-', 1, 8));
                obj.lowFrequencyModeShapes(Modal);
            end
            if params.animation
                fprintf('%s Mode shapes 3D animation\n', repmat('-', 1, 8));
                obj.lowFrequencyModeShapeMovie(Modal);
            end
            if params.Usteady
                fprintf('%s Steady-state damped response\n', repmat('-', 1, 8));
                obj.modeResponse(Modal.steady, "Modal", "Steady-State", "Displacement");
                obj.modeResponse(Modal.steady, "Modal", "Steady-State", "Velocity");
                obj.modeResponse(Modal.steady, "Modal", "Steady-State", "Acceleration");
                obj.modeResponse(Modal.steady, "Modal", "Steady-State", "Reaction");
            end
            if params.Utransient
                fprintf('%s Transient damped response\n', repmat('-', 1, 8));
                obj.modeResponse(Modal.transient, "Modal", "Transient", "Displacement");
            end
            if params.Uforced
                fprintf('%s Total forced damped response\n', repmat('-', 1, 8));
                obj.modeResponse(Modal.forced, "Modal", "Forced", "Displacement");
            end
            if params.Ufree
                fprintf('%s Free damped response\n', repmat('-', 1, 8));
                obj.modeResponse(Modal.free, "Modal", "Free", "Displacement");
            end
        end

        function impedanceMatrixMethod(obj, params, Impedance, Modal)
            fprintf('%s Plotting impedance results\n', repmat('-', 1, 6));
            obj.f = Impedance.w / 2 / pi;
            obj.modes2plot = params.modes2plot;
            obj.fn_harmonic = Impedance.fn_harmonic;
            if params.Usteady
                fprintf('%s Steady-state response\n', repmat('-', 1, 8));
                obj.modeResponse(Impedance.steady, "Impedance", "Steady-State", "Displacement");
            end
            if nargin == 4 && params.versus_modal
                fprintf('%s Impedance vs. Modal response\n', repmat('-', 1, 8));
                obj.compareResults(Impedance.steady, Modal.steady, "Displacement");
            end
        end

        function craigBamptonMethod(obj, params, CB)
            fprintf('%s Plotting Craig-Bampton results\n', repmat('-', 1, 6));
            if params.frequencies
                obj.frequencyDistribution(CB, obj.folder('Craig_Bampton'));
            end
        end
    end

    methods (Access = private)
        function finiteElementsMesh(obj, FEM)
            fig = figure('Units', 'inches', 'Position', [1, 1, 3.5, 2.8], 'Visible', 'off');
            axis equal;
            hold on;
            % Loop through each QUAD element and plot it
            for i = 1:size(FEM.mesh.connectivity, 1)
                elementNodes = FEM.mesh.connectivity(i, :);
                x = FEM.mesh.nodes(elementNodes, 1);
                y = FEM.mesh.nodes(elementNodes, 2);
                z = FEM.mesh.nodes(elementNodes, 3);
                patch(x, y, z, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'b', 'LineWidth', 1.5);
            end
            hold off;
            grid on;          
            xlabel('X (mm)');
            ylabel('Y (mm)');
            zlabel('Z (mm)');
            view(3);
            legend(sprintf('%d NASTRAN CQUAD4 Elements', FEM.mesh.elements.Ne), 'Box', 'off', 'Position', [0.45, 0.8, 0.1, 0.1]);
            exportgraphics(fig, fullfile(obj.folder('FEM'), 'FE_mesh.jpg'),'Resolution', 300);
            close(fig);
        end

        function finiteElementsMeshNodes(obj, FEM)
            fig = figure('Units', 'inches', 'Position', [1, 1, 3.5, 2.8], 'Visible', 'off');
            axis equal;
            hold on;
            scatter3(FEM.X(FEM.nng), FEM.Y(FEM.nng), FEM.Z(FEM.nng), 5, '*', 'red', 'LineWidth', 0.25);
            scatter3(FEM.X(FEM.nnf), FEM.Y(FEM.nnf), FEM.Z(FEM.nnf), 5, '*', 'blue', 'LineWidth', 0.25);
            quiver3(FEM.mesh.force.r(1), FEM.mesh.force.r(2), FEM.mesh.force.r(3), ...
                    FEM.mesh.force.F(1), FEM.mesh.force.F(2), FEM.mesh.force.F(3), 'MaxHeadSize', 10, 'LineWidth', 0.5, 'Color', 'k');
            hold off;
            grid on;          
            xlabel('X (mm)');
            ylabel('Y (mm)');
            zlabel('Z (mm)');
            view(3);
            legend({'Ground nodes', 'Free nodes', 'Force'}, 'Box', 'off', 'Position', [0.75, 0.8, 0.1, 0.1]);
            exportgraphics(fig, fullfile(obj.folder('FEM'), 'FE_mesh_nodes.jpg'),'Resolution', 300);
            close(fig);
        end

        function frequencyDistribution(obj, Model, path)
            % Plot the density function of natural frequencies
            % Ignore infinity values and then compute the histogram
            freqs = Model.fn(~isinf(Model.fn));
            freqs = freqs(freqs>0);
            [counts, bin_edges] = histcounts(freqs, 100, 'Normalization', 'pdf');
            bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
            % Estimate the distribution using a kernel density estimate (KDE)
            kernels = ["normal", "box", "triangle", "epanechnikov"];
            [kde, xi] = ksdensity(freqs, 'Support', [0, max(bin_edges)], 'Function','pdf', 'Kernel', kernels(1));

            % Plot the PDF histogram and KDE
            fig = figure('Units', 'inches', 'Position', [1, 1, 3.5, 2.8], 'Visible', 'off');
            bar(bin_centers, counts, 'hist');
            hold on;
            plot(xi, kde, 'r-', 'LineWidth', 2);
            xlabel('$f_n$ (Hz)');
            ylabel('$P(f_n)$');
            legend({'$f_n$', 'Kernel Density Estimate'}, 'Box', 'off');
            exportgraphics(fig, fullfile(path, 'P_fn.jpg'),'Resolution', 300);
            close(fig);
        end

        function lowFrequencyModeShapes(obj, Modal)
            % Function to plot the mode shapes and store them in a specified folder
            for mode = 1:length(Modal.zf(1,:))
                fig = figure('Units', 'inches', 'Position', [1, 1, 3.5, 3.5], 'Visible', 'off');
                %axis equal;
                hold on;
                scatter3(Modal.xg, Modal.yg, Modal.zg, '*', 'red', 'LineWidth', 0.5);
                scatter3(Modal.xf, Modal.yf, Modal.zf(:,mode), '*','blue', 'LineWidth', 0.5);
                grid on;
                title(sprintf('Mode shape %d', mode));
                xlabel('X (mm)');
                ylabel('Y (mm)');
                zlabel('Mode Shape');
                view(3);
                legend({'Ground nodes', 'Free nodes'}, 'Box', 'off', 'Position', [0.75, 0.8, 0.1, 0.1]);
                figname = sprintf('Mode_Shape_%d.jpg', mode);
                exportgraphics(fig, fullfile(obj.folder('Modal'), figname),'Resolution', 300);
                close(fig);
            end
        end

        function lowFrequencyModeShapeMovie(obj, Modal)
            % Create a 3D animation of the mode vibrations
            for mode = 1:length(Modal.zf(1,:))
                movie_name = sprintf('3D_Animation_Mode_%d.mp4', mode);
                movie_path = fullfile(obj.folder('Modal'), movie_name);
                % Initialize a video writer object
                v = VideoWriter(movie_path, 'MPEG-4');
                v.FrameRate = 24; v.Quality = 100;
                open(v);
                fig = figure('Units', 'inches', 'Position', [1, 1, 3.5, 3.5], 'Visible', 'off');
                %axis equal;
                % Generate and save frames
                wT = linspace(0, 2*pi, 96);

                for t = 1:length(wT)
                    clf;
                    hold on;
                    scatter3(Modal.xg, Modal.yg, Modal.zg, '*','red', 'LineWidth', 0.5);
                    scatter3(Modal.xf, Modal.yf, Modal.zf(:,mode)*cos(wT(t)), '*','blue', 'LineWidth', 0.5);
                    grid on;
                    colormap(jet);
                    zlim([-1, 1]);
                    xlabel('X (mm)');
                    ylabel('Y (mm)');
                    zlabel('Mode Shape');
                    title(sprintf('Mode Shape %d', mode));
                    view(3);
                    frame = getframe(fig);
                    writeVideo(v, frame);
                    pause(0.05);
                end
            end
            close(fig);
            close(v);
        end

        function modeResponse(obj, Model, Method, Mode, Type)
            [dt, nodes, UVAR] = obj.getSpecifiedResponse(Model, Type);       
            for node = 1:length(nodes)
                % Response for the node
                Resp = UVAR(node);
                % Node number and its DoF
                nn_dof = split(nodes(node), '-');            
                % Plot displacement, amplitude, and phase responses
                obj.fullResponse(Method, Mode, Type, nn_dof, dt, Resp.x);
                if Mode == "Steady-State"
                    obj.amplitudeResponse(Method, Mode, Type, nn_dof, Resp.A);
                    obj.phaseResponse(Method, Mode, Type, nn_dof, Resp.q);
                end
            end
        end

        function [dt, nodes, Response] = getSpecifiedResponse(obj, Model, Type)
            dt    = Model.delta_t;
            nodes = keys(Model.Ut);
            if Type == "Displacement"
                Response = values(Model.Ut);
            elseif Type == "Velocity"
                Response = values(Model.Vt);
            elseif Type == "Acceleration"
                Response = values(Model.At);
            else
                Response = values(Model.Rt);
            end
        end

        function fullResponse(obj, Method, Mode, Type, nn_dof, dt, Ut)
            fig = figure('Units', 'inches', 'Position', [1, 1, 3.5, 2.8], 'Visible', 'off');
            for mode = obj.modes2plot
                plot(dt, Ut(:, mode), 'DisplayName', sprintf('Mode %d', mode));
                hold on;
            end
            grid on;
            title(sprintf('Node %s DoF %s', nn_dof(1), nn_dof(2)));
            xlabel('Time (s)');
            ylabel(sprintf('%s Response %s', Type, obj.units(Type)));
            legend('Location','northeast', 'Box', 'on');

            figname = sprintf('%s_%s_Response_Node_%s_DoF_%s.jpg', Mode, Type, nn_dof(1), nn_dof(2));
            exportgraphics(fig, fullfile(obj.folder(Method), figname), 'Resolution', 300);
            close(fig);
        end

        function amplitudeResponse(obj, Method, Mode, Type, nn_dof, Amp)
            fig = figure('Units', 'inches', 'Position', [1, 1, 3.5, 2.8], 'Visible', 'off');
            semilogy(obj.f, Amp, 'DisplayName', 'Amplitude');
            hold on;
            y = min(nonzeros(Amp)) * ones(length(obj.fn_harmonic), 1);
            scatter(obj.fn_harmonic, y, 'g','^', 'DisplayName', 'Natural Frequency');
            hold off;
            grid on;
            title(sprintf('Node %s DoF %s', nn_dof(1), nn_dof(2)));
            xlabel('Frequency (Hz)');
            ylabel(sprintf('%s Amplitude %s', Type, obj.units(Type)));
            legend('Location','best', 'Box', 'off');

            figname = sprintf('%s_%s_Amplitude_Node_%s_DoF_%s.jpg', Mode, Type, nn_dof(1), nn_dof(2));
            exportgraphics(fig, fullfile(obj.folder(Method), figname), 'Resolution', 300);
            close(fig);
        end

        function phaseResponse(obj, Method, Mode, Type, nn_dof, Phase)
            fig = figure('Units', 'inches', 'Position', [1, 1, 3.5, 2.8], 'Visible', 'off');
            plot(obj.f, Phase, 'DisplayName','Phase');
            hold on;
            y = min(Phase) * ones(length(obj.fn_harmonic), 1);
            scatter(obj.fn_harmonic, y, 'g','^', 'DisplayName','Natural Frequency');
            hold off;
            grid on;
            title(sprintf('Node %s DoF %s', nn_dof(1), nn_dof(2)));
            xlabel('Frequency (Hz)');
            ylim([0, pi]);
            yticks([0, pi/6, pi/3, pi/2, 2*pi/3, 5*pi/6, pi]);
            yticklabels({'$0$', '$\pi/6$', '$\pi/3$', '$\pi/2$', '$2\pi/3$', '$5\pi/6$', '$\pi$'});
            ylabel(sprintf('%s Phase (rad)', Type));
            legend('Location','best', 'Box', 'off');

            figname = sprintf('%s_%s_Phase_Node_%s_DoF_%s.jpg', Mode, Type, nn_dof(1), nn_dof(2));
            exportgraphics(fig, fullfile(obj.folder(Method), figname), 'Resolution', 300);
            close(fig);
        end

        function compareResults(obj, Impedance, Modal, Type)
            nodes = keys(Modal.Ut);
            modal = values(Modal.Ut);
            impedance = values(Impedance.Ut);

            for node = 1:length(nodes)
                % Node number and its DoF
                nn_dof = split(nodes(node), '-');
                fig = figure('Name', 'Comparison', 'Units', 'inches', 'Position', [1, 1, 3.5, 5.6], 'Visible', 'off');

                subplot(2, 1, 1);
                semilogy(obj.f, modal(node).A, 'DisplayName', 'Modal');
                hold on;
                semilogy(obj.f, impedance(node).A, 'DisplayName', 'Impedance');
                hold on;
                % Display natural frequencies at the bottom of the plot
                y = min(modal(node).A) * ones(length(obj.fn_harmonic), 1);
                scatter(obj.fn_harmonic, y, 'g','^', 'DisplayName', 'Natural Frequency');
                hold off;
                grid on;
                title(sprintf('Node %s DoF %s', nn_dof(1), nn_dof(2)));
                xlabel('Frequency (Hz)');
                ylabel(sprintf('%s Amplitude %s', Type, obj.units(Type)));
                legend('Location', 'best', 'Box', 'off');

                subplot(2, 1, 2);
                plot(obj.f, modal(node).q, 'DisplayName', 'Modal');
                hold on;
                plot(obj.f, impedance(node).q, 'DisplayName', 'Impedance');
                hold on;
                y = min(modal(node).q) * ones(length(obj.fn_harmonic), 1);
                scatter(obj.fn_harmonic, y, 'g','^', 'DisplayName', 'Natural Frequency');
                hold off;
                grid on;
                xlabel('Frequency (Hz)');
                ylim([0, pi]);
                yticks([0, pi/6, pi/3, pi/2, 2*pi/3, 5*pi/6, pi]);
                yticklabels({'$0$', '$\pi/6$', '$\pi/3$', '$\pi/2$', '$2\pi/3$', '$5\pi/6$', '$\pi$'});
                ylabel(sprintf('%s Phase (rad)', Type));
                legend('Location', 'best', 'Box', 'off');

                figname = sprintf('Modal_vs_Impedance_Node_%s_DoF_%s.jpg', nn_dof(1), nn_dof(2));
                exportgraphics(fig, fullfile('Results/', figname), 'Resolution', 300);
                close(fig);
            end
        end
    end
end
