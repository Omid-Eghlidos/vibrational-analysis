classdef Plotting
    properties (Access = private)
        folder = dictionary('Modal', fullfile(pwd, 'Results/Modal'), ...
                            'Impedance', fullfile(pwd, 'Results/Impedance'));
        % Plotting modes
        plotting_modes = [];
        % All the excitation (w, f) and harmonic (fn_harmonic) frequencies
        w = []; f = []; fn_harmonic = [];
    end

    methods (Access = public)
        function [obj] = Plotting(plottin_modes, Modal)
            disp('---- Plotting');
            disp('------ Making folders for results');
            vals = values(obj.folder);
            for i = 1:length(keys(obj.folder))
                mkdir(vals(i));
            end
            obj.plotting_modes = plottin_modes;
            obj.fn_harmonic = Modal.fn_harmonic;
            obj.w = Modal.w;
            obj.f = Modal.w / 2 / pi;
        end

        function modalExpansion(obj, Modal, animation)
            disp('-- Plotting modal...')
            obj.frequencyDistribution(Modal);
            obj.lowFrequencyModeShapes(Modal);
            if nargin > 2
                obj.lowFrequencyModeShapeMovie(Modal);
            end
        end

        function dampedSteadyStateResponse(obj, Modal)
            disp('-- Plotting steady-state...');
            Model = Modal.steady;
            obj.modeResponse(Model, "Modal", "Steady-State", "Displacement");
            obj.modeResponse(Model, "Modal", "Steady-State", "Velocity");
            obj.modeResponse(Model, "Modal", "Steady-State", "Acceleration");
            obj.modeResponse(Model, "Modal", "Steady-State", "Reaction");
        end

        function dampedTransientResponse(obj, Modal)
            disp('-- Plotting transient...');
            Model = Modal.transient;
            obj.modeResponse(Model, "Modal", "Transient", "Displacement")
        end

        function dampedForcedResponse(obj, Modal)
            disp('-- Plotting forced...')
            Model = Modal.forced;
            obj.modeResponse(Model, "Modal", "Forced", "Displacement")
        end

        function dampedFreeResponse(obj, Modal)
            disp('-- Plotting free...');
            Model = Modal.free;
            obj.modeResponse(Model, "Modal", "Free", "Displacement");
        end

        function impedanceMatrixResponse(obj, Impedance)
            disp('-- Plotting impedance matrix...');
            Model = Impedance.steady;
            obj.modeResponse(Model, "Impedance", "Steady-State", "Displacement");
        end

        function compareModalAndImpedanceResults(obj, Modal, Impedance)
            disp('-- Plotting comparison of modal and impedance results');
            obj.compareResults(Modal.steady, Impedance.steady);
        end
    end

    methods (Access = private)
        function frequencyDistribution(obj, Modal)
            % Plot the density function of natural frequencies
            disp('---- Frequency distribution');
            % Ignore infinity values and then compute the histogram
            freqs = Modal.fn(~isinf(Modal.fn));
            [counts, bin_edges] = histcounts(freqs, 100, 'Normalization', 'pdf');
            bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
            % Estimate the distribution using a kernel density estimate (KDE)
            kernels = ["normal", "box", "triangle", "epanechnikov"];
            [kde, xi] = ksdensity(freqs, 'Support', [0, max(bin_edges)], 'Function','pdf', 'Kernel', kernels(1));

            % Plot the PDF histogram and KDE
            fig = figure(1);
            fig.Visible = 'off';
            bar(bin_centers, counts, 'hist');
            hold on;
            plot(xi, kde, 'r-', 'LineWidth', 2);
            xlabel('$f_n$ (Hz)');
            ylabel('$P(f_n)$');
            title('Distribution of Natural Frequencies');
            legend('$f_n$', 'Kernel Density Estimate');
            exportgraphics(fig, fullfile(obj.folder('Modal'), 'P_fn.jpg'),'Resolution', 300);
            close(fig);
        end

        function lowFrequencyModeShapes(obj, Modal)
            % Function to plot the mode shapes and store them in a specified folder
            disp('---- Mode shapes');
            for mode = 1:length(Modal.zf(1,:))
                fig = figure(mode);
                fig.Visible = "off";
                scatter3(Modal.xg, Modal.yg, Modal.zg, '*','red', 'LineWidth',0.5);
                hold on;
                scatter3(Modal.xf, Modal.yf, Modal.zf(:,mode), '*','blue', 'LineWidth',0.5);
                grid on;
                title(sprintf('Mode shape %d', mode));
                xlabel('x');
                ylabel('y');
                zlabel('Mode Shape');
                lgd = legend('Ground nodes', 'Free nodes');
                lgd.Position = [0.75, 0.8, 0.1, 0.1];
                legend show;
                figname = sprintf('Mode_Shape_%d.jpg', mode);
                exportgraphics(fig, fullfile(obj.folder('Modal'), figname),'Resolution', 300);
                close(fig);
            end
        end

        function lowFrequencyModeShapeMovie(obj, Modal)
            % Create a 3D animation of the mode vibrations
            disp('---- Creating 3D animation of mode shapes...');
            for mode = 1:length(Modal.zf(1,:))
                movie_name = sprintf('3D_Animation_Mode_%d.mp4', mode);
                movie_path = fullfile(obj.folder('Modal'), movie_name);
                % Initialize a video writer object
                v = VideoWriter(movie_path, 'MPEG-4');
                v.FrameRate = 24;
                open(v);
                fig = figure(mode);
                fig.Visible = "off";
                % Generate and save frames
                wT = linspace(0, 2*pi, 96);
                for t = 1:length(wT)
                    clf;
                    scatter3(Modal.xg, Modal.yg, Modal.zg, '*','red', 'LineWidth',0.5);
                    hold on;
                    scatter3(Modal.xf, Modal.yf, Modal.zf(:,mode)*cos(wT(t)), '*','blue', 'LineWidth',0.5);
                    grid on;
                    colormap(jet);
                    zlim([-1, 1]);
                    xlabel('X');
                    ylabel('Y');
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
            fprintf('---- %s\n', Type);
            [dt, nodes, UVAR] = obj.getSpecifiedResponse(Model, Type);

            for node = 1:length(nodes)
                % Response for the node
                Resp = UVAR(node);
                % Node number and its DoF
                nn_dof = split(nodes(node), '-');
                % Determine position, size, and number of figures
                position = [1, 1, 3.8, 2.5];
                num_figures = 1;
                if Mode == "Steady-State"
                    position = [1, 1, 3.8, 7.5];
                    num_figures = 3;
                end

                fig = figure('Name', Type, 'Units', 'inches', 'Position', position, "Visible", "off");
                subplot(num_figures, 1, 1);
                legends = [];
                for mode = obj.plotting_modes
                    plot(dt, Resp.x(:, mode))
                    legends = [legends, sprintf("Mode %d", mode)];
                    hold on;
                end
                grid on;
                title(sprintf('Node %s Degree of Freedom %s', nn_dof(1), nn_dof(2)));
                xlabel('Time (s)');
                ylabel(Type);
                legend(legends);

                if Mode == "Steady-State"
                    subplot(num_figures, 1, 2);
                    semilogy(obj.f, Resp.A);
                    hold on;
                    scatter(obj.fn_harmonic, min(Resp.A) * ones(length(obj.fn_harmonic), 1), 'g','^');
                    hold off;
                    grid on;
                    xlabel('Frequency (Hz)');
                    ylabel('Amplitude');
                    legend('Amplitude', 'Natural Frequency');
                    legend show;

                    subplot(3, 1, 3);
                    plot(obj.f, Resp.q);
                    hold on;
                    scatter(obj.fn_harmonic, min(Resp.q) * ones(length(obj.fn_harmonic), 1), 'g','^');
                    hold off;
                    grid on;
                    xlabel('Frequency (Hz)');
                    ylim([0, pi]);
                    yticks([0, pi/6, pi/3, pi/2, 2*pi/3, 5*pi/6, pi]);
                    yticklabels({'$0$', '$\pi/6$', '$\pi/3$', '$\pi/2$', '$2\pi/3$', '$5\pi/6$', '$\pi$'});
                    ylabel('Phase (rad)');
                    legend('Phase', 'Natural Frequency');
                    legend show;
                end

                figname = sprintf('%s_%s_Node_%s_DoF_%s.jpg', Mode, Type, nn_dof(1), nn_dof(2));
                exportgraphics(fig, fullfile(obj.folder(Method), figname), 'Resolution', 300);
                close(fig);
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

        function compareResults(obj, Modal, Impedance)
            nodes = keys(Modal.Ut);
            modal = values(Modal.Ut);
            impedance = values(Impedance.Ut);

            for node = 1:length(nodes)
                % Node number and its DoF
                nn_dof = split(nodes(node), '-');
                fig = figure('Name', 'Comparison', 'Units', 'inches', 'Position', [1, 1, 3.8, 5.0], "Visible", "off");

                subplot(2, 1, 1);
                semilogy(obj.f, modal(node).A);
                hold on;
                semilogy(obj.f, impedance(node).A);
                hold on;
                scatter(obj.fn_harmonic, min(modal(node).A) * ones(length(obj.fn_harmonic), 1), 'g','^');
                hold off;
                grid on;
                xlabel('Frequency (Hz)');
                ylabel('Amplitude');
                legend('Modal Expansion', 'Impedance Matrix', 'Natural Frequency');
                legend show;

                subplot(2, 1, 2);
                plot(obj.f, modal(node).q);
                hold on;
                plot(obj.f, impedance(node).q);
                hold on;
                scatter(obj.fn_harmonic, min(modal(node).q) * ones(length(obj.fn_harmonic), 1), 'g','^');
                hold off;
                grid on;
                xlabel('Frequency (Hz)');
                ylim([0, pi]);
                yticks([0, pi/6, pi/3, pi/2, 2*pi/3, 5*pi/6, pi]);
                yticklabels({'$0$', '$\pi/6$', '$\pi/3$', '$\pi/2$', '$2\pi/3$', '$5\pi/6$', '$\pi$'});
                ylabel('Phase (rad)');
                legend('Modal Expansion', 'Impedance Matrix', 'Natural Frequency');
                legend show;

                figname = sprintf('Modal_vs_Impedance_Node_%s_DoF_%s.jpg', nn_dof(1), nn_dof(2));
                exportgraphics(fig, fullfile('Results/', figname), 'Resolution', 300);
                close(fig);
            end
        end
    end
end
