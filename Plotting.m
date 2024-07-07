classdef Plotting
    properties (Access = private)
        folder = dictionary('Modal', fullfile(pwd, 'modal_results'), ...
            'Impedance', fullfile(pwd, 'impedance_results'));
    end

    methods (Access = public)
        function [obj] = Plotting()
            disp('---- Plotting');
            vals = values(obj.folder);
            for i = 1:length(keys(obj.folder))
                mkdir(vals(i));
            end
            disp('------ Made folder for results of each method');
        end

        function modalExpansion(obj, Modal)
            disp('-- Plotting modal properties...')
            obj.frequencyDistribution(Modal);
            %obj.lowFrequencyModeShapes(Modal);
            %obj.lowFrequencyModeShapeMovie(Modal);
        end

        function dampedSteadyStateResponse(obj, Modal)
            disp('-- Plotting steady-state response...')
            obj.modeResponse(Modal, "Steady-State", "Displacement")
            obj.modeResponse(Modal, "Steady-State", "Velocity");
            obj.modeResponse(Modal, "Steady-State", "Acceleration");
            obj.modeResponse(Modal, "Steady-State", "Reaction");
        end

        function dampedTransientResponse(obj, Modal)
            disp('-- Plotting transient response...')
            obj.modeResponse(Modal, "Transient", "Displacement")
            obj.modeResponse(Modal, "Transient", "Velocity");
            %obj.modeResponse(Modal, "Transient", "Acceleration");
        end

        function dampedForcedResponse(obj, Modal)
            disp('-- Plotting total forced response...')
            obj.modeResponse(Modal, "Forced", "Displacement")
        end

        function dampedFreeResponse(obj, Modal)
            disp('-- Plotting free response...')
            obj.modeResponse(Modal, "Free", "Displacement")
            obj.modeResponse(Modal, "Free", "Velocity");
        end

        function impedanceMatrixResponse(obj, Impedance)
            disp('-- Plotting impedance matrix response...')
            obj.modeResponse(Impedance, "Impedance", "Displacement")
        end
    end

    methods (Access = private)
        function frequencyDistribution(obj, Modal)
            % Plot the density function of natural frequencies
            disp('---- Frequency distribution');
            % Ignore infinity values
            freqs = Modal.fn(~isinf(Modal.fn));
            % Compute the histogram
            [counts, bin_edges] = histcounts(freqs, 100, 'Normalization', 'pdf');
            % Calculate bin centers
            bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

            % Plot the PDF
            fig = figure(1);
            bar(bin_centers, counts, 'hist');
            hold on;
            % Plot a smooth version using a kernel density estimate
            kernels = ["normal", "box", "triangle", "epanechnikov"];
            [fn, xi] = ksdensity(freqs, 'Support', [0, max(bin_edges)], 'Function','pdf', 'Kernel', kernels(1));
            plot(xi, fn, 'r-', 'LineWidth', 2);
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
                scatter3(Modal.xg, Modal.yg, Modal.zg, '*','red', 'LineWidth',0.5);
                hold on;
                scatter3(Modal.xf, Modal.yf, Modal.zf(:,mode), '*','blue', 'LineWidth',0.5);
                grid on;
                title(['Mode shape ', num2str(mode)])
                xlabel('x')
                ylabel('y')
                zlabel('Mode Shape')
                lgd = legend('Ground nodes', 'Free nodes');
                lgd.Position = [0.75, 0.8, 0.1, 0.1];
                legend show
                figname = sprintf('Mode_Shape_%d.jpg', mode);
                exportgraphics(figure(mode), fullfile(obj.folder('Modal'), figname),'Resolution', 300);
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
                % Create a figure for the animation
                fig = figure(1);
                % Generate and save frames
                wT = linspace(0, 2*pi, 96);
                for t = 1:length(wT)
                    clf;
                    % Update the surface data
                    scatter3(Modal.xg, Modal.yg, Modal.zg, '*','red', 'LineWidth',0.5);
                    hold on;
                    scatter3(Modal.xf, Modal.yf, Modal.zf(:,mode)*cos(wT(t)), '*','blue', 'LineWidth',0.5);
                    grid on;
                    colormap(jet);
                    %colorbar;
                    zlim([-1, 1]);
                    xlabel('X');
                    ylabel('Y');
                    zlabel('Mode Shape');
                    title(['Mode Shape ', num2str(mode)]);
                    view(3);
                    frame = getframe(fig);
                    writeVideo(v, frame);
                    pause(0.05);
                end
            end
            close(fig);
            close(v);
        end

        function modeResponse(obj, Model, Mode, Type)
            fprintf('---- %s %s\n', Mode, Type);
            if Mode ~= "Forced"
                [dt, nodes, UVAR] = obj.getSpecifiedResponse(Model, Mode, Type);
            else
                [dt, nodes, Up] = obj.getSpecifiedResponse(Model, "Transient", Type);
                [ ~,     ~, Uh] = obj.getSpecifiedResponse(Model, "Steady-State", Type);
            end

            for node = 1:length(nodes)
                if Mode ~= "Forced"
                    Resp = UVAR(node);
                    modes = [2, 3, 15];
                else
                    Resp.x = Up(node).x + Uh(node).x;
                    Resp.A = Up(node).A + Uh(node).A;
                    Resp.q = Uh(node).q;
                    modes = [2, 15];
                end

                % Node number and its DoF
                nn_dof = split(nodes(node), '-');

                fig = figure('Name', Type, 'Units', 'inches', 'Position', [1, 1, 3.8, 7.5]);
                subplot(3, 1, 1);
                legends = [];
                for mode = modes
                    plot(dt, Resp.x(:, mode))
                    if mode == 2
                        legends = [legends, sprintf("Mode = %d", (mode-1))];
                    else
                        legends = [legends, sprintf("Mode = %d", (mode))];
                    end
                    hold on;
                end
                grid on;
                title(sprintf('Node %s Degree of Freedom %s', nn_dof(1), nn_dof(2)));
                xlabel('Time (s)');
                ylabel(Type);
                legend(legends);

                subplot(3, 1, 2)
                f = Model.w / 2 / pi;
                if Mode == "Free"
                    f = Model.fn_free;
                end
                semilogy(f, Resp.A)
                hold on;
                scatter(Model.fn_harmonic, min(Resp.A) * ones(length(Model.fn_harmonic), 1), 'g','^');
                hold off;
                grid on;
                xlabel('Frequency (Hz)');
                ylabel('Amplitude');
                legend('Amplitude', 'Natural Frequency')
                legend show

                subplot(3, 1, 3)
                plot(f, Resp.q)
                hold on;
                scatter(Model.fn_harmonic, min(Resp.q) * ones(length(Model.fn_harmonic), 1), 'g','^');
                hold off;
                grid on;
                xlabel('Frequency (Hz)');
                ylim([0, pi]);
                yticks([0, pi/6, pi/3, pi/2, 2*pi/3, 5*pi/6, pi]);
                yticklabels({'$0$', '$\pi/6$', '$\pi/3$', '$\pi/2$', '$2\pi/3$', '$5\pi/6$', '$\pi$'});
                ylabel('Phase (rad)');
                legend('Phase', 'Natural Frequency')
                legend show

                figname = sprintf('%s_%s_Node_%s_DoF_%s.jpg', Mode, Type, nn_dof(1), nn_dof(2));
                exportgraphics(fig, fullfile(obj.folder('Modal'), figname), 'Resolution', 300);
                close(fig);
            end
        end

        function [dt, nodes, Response] = getSpecifiedResponse(obj, Model, Mode, Type)
            if Mode == "Steady-State"
                dt    = Model.steady.delta_t;
                nodes = keys(Model.steady.Ut);
                if Type == "Displacement"
                    Response = values(Model.steady.Ut);
                elseif Type == "Velocity"
                    Response = values(Model.steady.Vt);
                elseif Type == "Acceleration"
                    Response = values(Model.steady.At);
                else
                    Response = values(Model.steady.Rt);
                end
            elseif Mode == "Transient"
                dt    = Model.transient.delta_t;
                nodes = keys(Model.transient.Ut);
                if Type == "Displacement"
                    Response = values(Model.transient.Ut);
                elseif Type == "Velocity"
                    Response = values(Model.transient.Vt);
                elseif Type == "Acceleration"
                    Response = values(Model.transient.At);
                end
            elseif Mode == "Free"
                dt    = Model.free.delta_t;
                nodes = keys(Model.free.Ut);
                if Type == "Displacement"
                    Response = values(Model.free.Ut);
                elseif Type == "Velocity"
                    Response = values(Model.free.Vt);
                end
            elseif Mode == "Impedance"
                dt    = Model.steady.delta_t;
                nodes = keys(Model.steady.Ut);
                if Type == "Displacement"
                    Response = values(Model.steady.Ut);
                end
            end
        end
    end
end
