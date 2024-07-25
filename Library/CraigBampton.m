classdef CraigBampton
    properties (Access = public)
        % Selected substructure
        substructures = dictionary();
    end

    properties (Access = private)
        % External boundary nodes of the whole structure
        external_boundary_nodes = [];
        % Keep track of selected elements for all the substrcutures
        selected_substructures_elements = [];
        % Craig-Bampton transformation matrix
        Phi_CB = [];
        % Craig-Bampton reduced mass and stiffness matrices
        Mcb = []; Kcb = [];
    end

    methods (Access = public)
        function [obj] = CraigBampton(FEM)
            % Find the external boundary nodes of the whole structure
            obj = obj.findExternalBoundaryNodes(FEM);
            % Keep track of selected elements in the whole structure
            obj.selected_substructures_elements = false(FEM.mesh.elements.Ne, 1);
        end

        function [obj] = addSubstructures(obj, FEM)
            % Open a figure of the model and ask the user to select the substructures
            add_substructure = true;
            substructure_id = 0;
            obj.plotStructuralMesh(FEM);
            while add_substructure
                substructure_id = substructure_id + 1;
                fprintf('%s Add substructure %d by drawing a rectangle\n', repmat('-', 1, 6), substructure_id);
                substructure = obj.selectSubstructureAndFindItsFeatures(FEM, substructure_id);
                obj.substructures(substructure_id) = substructure;
                obj.selected_substructures_elements(substructure.elements) = true;
                add_substructure = obj.askUserApproval("Add");
            end
            % If there are any unselected elements add them as the last substructure
            if any(~obj.selected_substructures_elements)
                substructure_id = substructure_id + 1;
                substructure = obj.findSubstructureFeatures(FEM, ~obj.selected_substructures_elements);
                obj.substructures(substructure_id) = substructure;
                obj.selected_substructures_elements(substructure.elements) = true;
            end
            fprintf('%s Continue analysis with a total of %d substructures\n', repmat('-', 1, 6), substructure_id);
            obj.plotStructuralMesh(FEM);
        end

        function [obj] = performCraigBamptonReduction(obj)
            % Perform Craig-Bampton reduction for selected substructures
        end
    end

    methods (Access = private)
        function [approved] = askUserApproval(obj, Type)
            % Ask user to approve the selected substructure or if to add another substructure
            if Type == "Verify"
                switch questdlg('Do you approve the selected substructure?', 'Verify Substructure', 'Yes', 'No', 'Yes')
                    case 'Yes'
                        approved = true;
                        fprintf('%s Store the selected substructure\n', repmat('-', 1, 8));
                    case 'No'
                        approved = false;
                        fprintf('%s Select another substructure\n', repmat('-', 1, 8));
                end
            elseif Type == "Add"
                if all(obj.selected_substructures_elements)
                    approved = false;
                else
                    substructure_id = max(keys(obj.substructures)) + 1;
                    switch questdlg('Do you want to add another substructure?', 'Add Substructure', 'Yes', 'No', 'Yes')
                        case 'Yes'
                            approved = true;
                            fprintf('%s Add substructure %d\n', repmat('-', 1, 8), substructure_id);
                        case 'No'
                            approved = false;
                    end
                end
            end
        end

        function [legends] = plotStructuralMesh(obj, FEM)
            % Plot structural mesh and if there are previously selected substructures
            close all;
            fig = figure();
            hold on;
            axis equal;
            % Ensure selected rectangle is on top of the mesh by offset z by structure thickness
            z_offset = FEM.mesh.dimension.z(2) - FEM.mesh.dimension.z(1);
            for i = 1:size(FEM.mesh.connectivity, 1)
                element_nodes = FEM.mesh.connectivity(i, :);
                x = FEM.mesh.nodes(element_nodes, 1);
                y = FEM.mesh.nodes(element_nodes, 2);
                z = FEM.mesh.nodes(element_nodes, 3) - z_offset;
                Structure = patch(x, y, z, [0 0.4470 0.7410], 'EdgeColor', 'k', 'LineWidth', 0.5);
            end
            grid on;
            xlims = xlim; ylims = ylim; zlims = zlim;
            xlim([xlims(1)-1.0, xlims(2)+1.0]);
            ylim([ylims(1)-1.0, ylims(2)+1.0]);
            zlim([zlims(1)-1.0, zlims(2)+1.0]);
            xlabel('X (mm)');
            ylabel('Y (mm)');
            zlabel('Z (mm)');
            view(2);

            % Add previously selected substructures if they exist
            legends = dictionary("Structure Mesh", Structure);
            if isConfigured(obj.substructures)
                for i = 1:length(keys(obj.substructures))
                    [Names, Plots, ~] = obj.plotSubstructure(FEM, obj.substructures(i), i);
                    legends(Names) = Plots;
                end
            end

            if all(obj.selected_substructures_elements)
                fig.Visible = 'off';
                Values = values(legends);
                [Keys, sortIdx] = sort(keys(legends), 'descend');
                legend(Values(sortIdx), num2cell(Keys));
                % Save figure of the substructures
                exportgraphics(gcf, fullfile(pwd, 'Results/Craig_Bampton/substructures.jpg'), 'Resolution', 300);
                close(fig);
            end
        end

        function [Names, Plots, Patches] = plotSubstructure(obj, FEM, substructure, substructure_id)
            % Plot the substructure elements and its interior and boundary nodes
            z_offset = FEM.mesh.dimension.z(2) - FEM.mesh.dimension.z(1);
            Patches = [];
            for i = 1:size(FEM.mesh.connectivity, 1)
                if substructure.elements(i)
                    element_nodes = FEM.mesh.connectivity(i, :);
                    x = FEM.mesh.nodes(element_nodes, 1);
                    y = FEM.mesh.nodes(element_nodes, 2);
                    z = FEM.mesh.nodes(element_nodes, 3) - z_offset;
                    Patch = patch(x, y, z, 'FaceColor', obj.substructureColor(substructure_id), 'EdgeColor', 'k', 'LineWidth', 0.5);
                    Patches = [Patches; Patch];
                end
            end

            Ni = substructure.nodes.I; Nb = substructure.nodes.B;
            interior = plot3(FEM.mesh.nodes(Ni, 1), FEM.mesh.nodes(Ni, 2), FEM.mesh.nodes(Ni, 3), 'go', 'MarkerSize', 8, 'LineWidth', 2);
            boundary = plot3(FEM.mesh.nodes(Nb, 1), FEM.mesh.nodes(Nb, 2), FEM.mesh.nodes(Nb, 3), 'ro', 'MarkerSize', 8, 'LineWidth', 2);

            % Return the legends for the substructure
            Names = [sprintf("Substructure %d", substructure_id), "Interior Nodes", "Boundary Nodes"];
            Plots = [Patches(end), interior, boundary];
        end

        function [substructure, legends] = selectSubstructureAndFindItsFeatures(obj, FEM, substructure_id)
            % Select a non-overlapping substructure by drawing a rectangle on the mesh and find its features
            approved = false;
            while ~approved
                overlap = true;
                while overlap
                    % Limit selection to the structure mesh area
                    drawing_area = [FEM.mesh.dimension.x(1), FEM.mesh.dimension.y(1), ...
                        diff(FEM.mesh.dimension.x), diff(FEM.mesh.dimension.y)];
                    rectangle = drawrectangle('Label', sprintf('Substructure %d', substructure_id), ...
                        'Color', obj.substructureColor(substructure_id), 'DrawingArea', drawing_area);
                    selected_nodes = find(inROI(rectangle, FEM.mesh.nodes(:,1), FEM.mesh.nodes(:,2)));
                    selected_elements_nodes = ismember(FEM.mesh.connectivity, selected_nodes);
                    selected_elements = false(size(FEM.mesh.connectivity, 1), 1);
                    for i = 1:size(FEM.mesh.connectivity, 1)
                        if any(selected_elements_nodes(i, :))
                            selected_elements(i) = true;
                        end
                    end
                    substructure = obj.findSubstructureFeatures(FEM, selected_elements);
                    overlap = obj.checkSubstructuresOverlap(substructure, rectangle);
                end
                [Names, Plots, Patches] = obj.plotSubstructure(FEM, substructure, substructure_id);
                legends = dictionary(Names, Plots);
                approved = obj.askUserApproval("Verify");
                if ~approved
                    delete(rectangle);
                    delete(Plots);
                    delete(Patches);
                end
            end
        end

        function [overlap] = checkSubstructuresOverlap(obj, substructure, rectangle)
            % Check if the selected substructure has any overlap with previously selected ones
            overlap = false;
            if isConfigured(obj.substructures)
                for i = 1:length(keys(obj.substructures))
                    if any(ismember(find(substructure.elements), find(obj.substructures(i).elements)))
                        fprintf('WARNING: Selected substructure overlap with another substructure(s), select again\n');
                        overlap = true;
                        delete(rectangle);
                        return;
                    end
                end
            end
        end

        function [substructure] = findSubstructureFeatures(obj, FEM, elements)
            % Find substructure's different nodes and its corresponding matrices
            % Find the substructure's different nodes
            [nodes, dofs] = obj.findSubstructureInteriorAndBoundaryNodes(FEM, elements);
            % Get the substructure's matrices from corresponding structural matrices
            [m, c, k, f] = obj.getSubstructureReducedMatrices(FEM, dofs);
            % Store substructure features
            substructure = struct('elements', elements , 'nodes', nodes, 'dofs', dofs, ...
                'm', m, 'c', c, 'k', k, 'f', f);
        end

        function [nodes, dofs] = findSubstructureInteriorAndBoundaryNodes(obj, FEM, elements)
            % Find the substructure's elements and interior and boundary nodes
            % Check neighboring elements to find boundary and interior nodes
            boundary = []; interior = [];
            for i = 1:size(FEM.mesh.connectivity, 1)
                if elements(i)
                    element_nodes = FEM.mesh.connectivity(i, :);
                    for j = 1:length(element_nodes)
                        % Check if the node is connected to any element outside of the ROI
                        neighboring_elements = find(any(FEM.mesh.connectivity == element_nodes(j), 2));
                        if any(~elements(neighboring_elements))
                            boundary = [boundary; element_nodes(j)];
                        else
                            interior = [interior; element_nodes(j)];
                        end
                    end
                end
            end
            % Remove any external boundary nodes in substructure's
            % interior nodes and add them to its boundary nodes
            interior = sort(unique(interior));
            external_nodes = intersect(interior, obj.external_boundary_nodes);
            boundary = sort(unique([unique(boundary); external_nodes]));
            interior(ismember(interior, external_nodes)) = [];
            nodes = struct('B', boundary, 'I', interior);

            % Boundary nodes DoFs and interior nodes DoFs
            dofs_boundary = reshape(repmat((boundary(:) - 1) * 6, 1, 6)' + (1:6)', [], 1);
            dofs_interior = reshape(repmat((interior(:) - 1) * 6, 1, 6)' + (1:6)', [], 1);
            dofs = struct('B', dofs_boundary, 'I', dofs_interior);
        end

        function [m, c, k, f] = getSubstructureReducedMatrices(obj, FEM, dofs)
            % Find the mass, damping, and stiffness matrices of the
            % substructure from the the corresponding structure's matrices
            % Substructure's mass block matrics
            m = struct();
            m.BB = FEM.M(dofs.B, dofs.B); m.BI = FEM.M(dofs.B, dofs.I);
            m.IB = FEM.M(dofs.I, dofs.B); m.II = FEM.M(dofs.I, dofs.I);
            % Substructure's damping block matrics
            c = struct();
            c.BB = FEM.C(dofs.B, dofs.B); c.BI = FEM.C(dofs.B, dofs.I);
            c.IB = FEM.C(dofs.I, dofs.B); c.II = FEM.C(dofs.I, dofs.I);
            % Substructure's stiffness block matrics
            k = struct();
            k.BB = FEM.K(dofs.B, dofs.B); k.BI = FEM.K(dofs.B, dofs.I);
            k.IB = FEM.K(dofs.I, dofs.B); k.II = FEM.K(dofs.I, dofs.I);
            % Substructure's force block vector
            f = struct();
            f.B = FEM.F(dofs.B);          f.I = FEM.F(dofs.I);
        end

        function [obj] = findExternalBoundaryNodes(obj, FEM)
            % Find the external boundary nodes of the entire structure
            % The coordinates of these nodes matches either x or y dimension
            % For QUAD elements they have only 3 neighbors
            node_count = zeros(size(FEM.mesh.nodes, 1), 1);
            for i = 1:size(FEM.mesh.connectivity, 1)
                node_count(FEM.mesh.connectivity(i, :)) = node_count(FEM.mesh.connectivity(i, :)) + 1;
            end
            obj.external_boundary_nodes = find(node_count < 4);
        end

        function [color] = substructureColor(obj, substructure_id)
            % Return a color depending on the id of the substructure
            %patch_colors = ["#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F"];
            colors = {[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], ...
                [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], ...
                [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840]};
            color = colors{substructure_id};
        end

        function [obj] = reducedMatrices(obj)
            % Craig-Bampton reduced mass and stiffness matrices
            obj.Mcb = obj.Phi_CB' * obj.M * obj.Phi_CB;
            obj.Kcb = obj.Phi_CB' * obj.K * obj.Phi_CB;
        end
    end
end
