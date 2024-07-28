classdef CraigBampton
    properties (Access = public)
        % Sorted eigenvalues and eigenvectors of the structure
        Phi = [], fn = [];
    end

    properties (Access = private)
        % Output file to save or load the substructures
        output_file = "Inputs/substructures.mat";
        % Setting for loading substructures if trur
        load_substructures = false;
        % Tolerance for assessing the results
        tolerance = 1e-3;
        % Number of modes to retain
        Nmodes = 10;
        % External boundary nodes of the whole structure
        external_boundary_nodes = [];
        % Keep track of selected elements for all the substrcutures
        selected_substructures_elements = [];
        % Selected substructure's physical and Craig-Bampton reduced properties
        substructures = dictionary();
        % Craig-Bampton whole structural system properties
        structure = struct('M', [], 'K', [], 'F', []);
    end

    methods (Access = public)
        function [obj] = CraigBampton(params)
            obj.load_substructures = params.use_existing_substructures;
            obj.output_file = params.output_file;
            obj.tolerance = params.tolerance;
            obj.Nmodes = params.Nmodes;
        end

        function [obj] = addSubstructures(obj, FEM)
            % Function to select substructures for Crain-Bampton analysis            
            if isfile(obj.output_file) && obj.load_substructures
                % Load substructures from "Inputs\substructures.mat" if exists 
                obj = obj.saveLoadSubstructures("Load");
                obj.selected_substructures_elements = true(FEM.mesh.elements.Ne, 1);
                obj.plotStructuralMesh(FEM);
            else
                % Select from the structural mesh plot
                obj = obj.findStructureExternalBoundaryNodes(FEM);
                obj.selected_substructures_elements = false(FEM.mesh.elements.Ne, 1);
                obj = obj.selectSubstructures(FEM);
            end
        end

        function [obj] = performCraigBamptonReduction(obj)
            % Perform Craig-Bampton reduction for selected substructures
            fprintf('%s Start Craig-Bampton reduction process\n', repmat('-', 1, 6));
            obj = obj.computeReducedSubstructureFeatures();
            obj = obj.assembleStructureGlobalMatrices();
            obj = obj.computeStructureModalProperties();
        end
    end

    methods (Access = private)
        function [obj] = selectSubstructures(obj, FEM)
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
                fprintf('%s Add the remaining elements as substructure %d\n', repmat('-', 1, 6), substructure_id);
                substructure = obj.findSubstructureFeatures(FEM, ~obj.selected_substructures_elements);
                obj.substructures(substructure_id) = substructure;
                obj.selected_substructures_elements(substructure.elements) = true;
            end
            fprintf('%s Continue analysis with a total of %d substructures\n', repmat('-', 1, 6), substructure_id);
            obj.plotStructuralMesh(FEM);
            % Save substructures into a file for later use
            obj.saveLoadSubstructures("Save");
        end

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

        function [obj] = saveLoadSubstructures(obj, Type)
            % Saves and loads substructures property from or into an output file
            if Type == "Save"
                fprintf('%s Saving substructures into %s\n', repmat('-', 1, 6), obj.output_file);
                % Convert substructures dictionary to struct and then save them into a file
                Keys = keys(obj.substructures);
                Values = values(obj.substructures);
                Substructures = struct();
                for i = 1:length(Keys)
                    Substructures.(sprintf('s%d', Keys(i))) = Values(i);
                end
                save(obj.output_file, 'Substructures');
            elseif Type == "Load"
                fprintf('%s Loading substructures from %s\n', repmat('-', 1, 6), obj.output_file);
                % Convert substructures dictionary to struct and then save them into a file
                Substructures = load(obj.output_file).Substructures;
                Fields = fieldnames(Substructures);                
                for i = 1:length(Fields)
                    obj.substructures(str2num(strip(Fields{i}, 's'))) = Substructures.(Fields{i});
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
                for i = 1:numEntries(obj.substructures)
                    [Names, Plots, ~] = obj.plotSubstructure(FEM, obj.substructures(i), i);
                    legends(Names) = Plots;
                end
            end

            if all(obj.selected_substructures_elements)
                fig.Visible = 'off';
                Values = values(legends);
                [Keys, sortIdx] = sort(keys(legends), 'descend');
                legend(Values(sortIdx), num2cell(Keys));
                % Save figure of the substructures if it is selected
                if ~obj.load_substructures
                    schematic_path = 'Results/Craig_Bampton/substructures.jpg';
                    exportgraphics(gcf, fullfile(pwd, schematic_path), 'Resolution', 300);
                end
                pause(3);
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
            % Plot different nodes on the mesh
            settings = struct('LineWidth', 2,  'MarkerSize', 8);
            NI = substructure.nodes.I;  NBd = substructure.nodes.Bd;  NBg = substructure.nodes.Bg;
            I = plot3(FEM.mesh.nodes(NI, 1), FEM.mesh.nodes(NI, 2), FEM.mesh.nodes(NI, 3), 'go', settings);
            Bd = plot3(FEM.mesh.nodes(NBd, 1), FEM.mesh.nodes(NBd, 2), FEM.mesh.nodes(NBd, 3), 'ro', settings);
            Bg = plot3(FEM.mesh.nodes(NBg, 1), FEM.mesh.nodes(NBg, 2), FEM.mesh.nodes(NBg, 3), 'ko', settings);

            % Return the legends for the substructure
            Names = [sprintf("Substructure %d", substructure_id), ...
                     "Interior Nodes", "Dependent Boundary Nodes", "Independent Boundary Nodes"];
            Plots = [Patches(end), I, Bd, Bg];
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
                for i = 1:numEntries(obj.substructures)
                    if any(ismember(find(substructure.elements), find(obj.substructures(i).elements)))
                        fprintf('WARNING: Overlap between substructures, select again\n');
                        overlap = true;
                        delete(rectangle);
                        return;
                    end
                end
            end
        end

        function [substructure] = findSubstructureFeatures(obj, FEM, elements)
            % Find substructure's different nodes and its corresponding matrices
            fprintf('%s Finding substructures physical features\n', repmat('-', 1, 8));
            nodes = obj.findSubstructureInteriorAndBoundaryNodes(FEM, elements);
            dofs = obj.findSubstructuresLocalAndGlobalDoFs(nodes);
            [m, k, F] = obj.getSubstructureMatrices(FEM, dofs);
            u = obj.computeInitialPhysicalResponse(dofs, k, F);
            % Store substructure's features
            substructure = struct('elements', elements , 'nodes', nodes, ...
                             'dofs', dofs, 'm', m, 'k', k, 'F', F, 'u', u);
        end

        function [nodes] = findSubstructureInteriorAndBoundaryNodes(obj, FEM, elements)
            % Find the substructure's elements and interior and boundary nodes
            % Check neighboring elements to find dependent boundary (Bd) and interior (I) nodes
            Bd = []; I = [];
            for i = 1:size(FEM.mesh.connectivity, 1)
                if elements(i)
                    element_nodes = FEM.mesh.connectivity(i, :);
                    for j = 1:length(element_nodes)
                        % Check if the node is connected to any element outside of the ROI
                        neighboring_elements = find(any(FEM.mesh.connectivity == element_nodes(j), 2));
                        if any(~elements(neighboring_elements))
                            Bd = [Bd; element_nodes(j)];
                        else
                            I = [I; element_nodes(j)];
                        end
                    end
                end
            end
            % Keep sorted unique nodes for each type of boundary and the interior
            I = sort(unique(I));      
            Bd = sort(unique(Bd));
            Bg = sort([intersect(Bd, obj.external_boundary_nodes); intersect(I, obj.external_boundary_nodes)]);
            Bd(ismember(Bd, Bg)) = [];
            % Remove any external/independent boundary from interior nodes   
            I(ismember(I, Bg)) = [];
            nodes = struct('Bd', Bd, 'Bg', Bg, 'B', unique([Bd; Bg]), 'I', I);
        end

        function [dofs] = findSubstructuresLocalAndGlobalDoFs(obj, nodes)
            % Find DoF number for for dependent (Bd), independednt (Bg),
            % and total boundary (B), and interior (I)
            % Global
            Bd = reshape(repmat((nodes.Bd(:) - 1) * 6, 1, 6)' + (1:6)', [], 1);
            Bg = reshape(repmat((nodes.Bg(:) - 1) * 6, 1, 6)' + (1:6)', [], 1);
            B  = reshape(repmat((nodes.B(:) - 1)  * 6, 1, 6)' + (1:6)', [], 1);
            I  = reshape(repmat((nodes.I(:) - 1)  * 6, 1, 6)' + (1:6)', [], 1);
            G = struct('Bd', Bd, 'Bg', Bg, 'B', B, 'I', I);

            % Local (relative within each type)
            [~, bd, ~] = intersect([B; I], Bd);
            [~, bg, ~] = intersect([B; I], Bg);
            [~,  b, ~] = intersect([B; I], B);
            [~,  i, ~] = intersect([B; I], I);
            L = struct('Bd', bd, 'Bg', bg, 'B', b, 'I', i);

            % Return both local and gloab dofs
            dofs = struct('g', G, 'l', L);
        end

        function [u] = computeInitialPhysicalResponse(obj, dofs, k, F)
            % Substructure total physical degrees of freedom
            U = [k.BB, k.BI; k.IB, k.II] \ [F.B; F.I];
            % Separated response for each type of DoF
            u = struct('Bd', U(dofs.l.Bd), 'Bg',  U(dofs.l.Bg), ...
                        'B', U(dofs.l.B),   'I',  U(dofs.l.I));
        end

        function [m, k, F] = getSubstructureMatrices(obj, FEM, dofs)
            % Find the mass, damping, and stiffness block matrices and force
            % block vector of the substructure from the the corresponding 
            % structure's matrices
            % Boundary and inetrior dofs to consider
            B = dofs.g.B;  I = dofs.g.I;
            % Mass block matrics
            m = struct('BB', FEM.M(B, B), 'BI', FEM.M(B, I),...
                       'IB', FEM.M(I, B), 'II', FEM.M(I, I));
            % Stiffness block matrics
            k = struct('BB', FEM.K(B, B), 'BI', FEM.K(B, I),...
                       'IB', FEM.K(I, B), 'II', FEM.K(I, I));
            % Force block vector
            F = struct('B', FEM.F(B), 'I', FEM.F(I));
        end

        function [obj] = findStructureExternalBoundaryNodes(obj, FEM)
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
            colors = {[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], ...
                [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], ...
                [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840]};
            color = colors{substructure_id};
        end

        function [obj] = computeReducedSubstructureFeatures(obj)
            % Craig-Bampton method for finding the substructure reduced matrices
            % Size of the globa matrices
            Ng = 0;
            for r = 1:numEntries(obj.substructures)
                fprintf('%s Computing generalized features for substructure %d ...\n', repmat('-', 1, 8), r);
                substructure = obj.substructures(r);
                phi_bar =  struct('C', [], 'N', []);
                [phi_bar.C] = obj.determineConstraintModes(substructure);
                % Use the first 15 modes, similar to modal expansion for comparison
                [w, phi_bar.N] = obj.computeSubstructuresGeneralizedNormalModes(substructure);
                alpha = obj.computeTransformationOfPhysicalToGeneralizedCoordinates(substructure, phi_bar);
                p = obj.findSubstructureGeneralizedCoordinates(substructure, alpha);
                [m_bar, k_bar, P] = obj.computeSubstructureGeneralizedMatrices(substructure, phi_bar, w, alpha);
                reduced_fields = {'phi_bar', 'w', 'p', 'alpha', 'm_bar', 'k_bar', 'P'};
                reduced_properties = {phi_bar, w, p, alpha, m_bar, k_bar, P};
                for i = 1:length(reduced_fields)
                    substructure.(reduced_fields{i}) = reduced_properties{i};
                end
                obj.substructures(r) = substructure;
                Ng = Ng + size(m_bar, 1);
            end
            % Initialize the structure (global) matrices
            obj.structure = struct('M', zeros(Ng), 'K', zeros(Ng), 'F', zeros(Ng, 1));
        end

        function [phi_bar_C] = determineConstraintModes(obj, substructure)
            % Find the constraint modes (mode shapes of the interior
            % freedoms) due to successive displacement of boundary points
            % Find the desired matrix of constraint      
            phi_bar_C = -substructure.k.II \ substructure.k.IB;
        end

        function [w_N, phi_bar_N] = computeSubstructuresGeneralizedNormalModes(obj, substructure)
            % Compute substructure normal modes for interior mass and 
            % stiffness matrices of the substructure in generalized coordinates p
            % Find substructure's eigenvalues and eigenvectors for interior nodes
            [phi_I, lambda] = eig(substructure.k.II, substructure.m.II);
            [w2, idx] = sort(diag(lambda));
            phi_I = phi_I(:, idx);
            w_I = sqrt(w2);
            % Frequencies and shapes for the specified modes
            w_N = w_I(1:obj.Nmodes);
            phi_bar_N = phi_I(:,1:obj.Nmodes);
        end

        function [alpha] = computeTransformationOfPhysicalToGeneralizedCoordinates(obj, substructure, phi_bar)
            % Find the alpha transformation between substructure's 
            % physical coordinates delta and its generalized coordinates p
            I = eye(size(substructure.m.BB,1), size(phi_bar.C,2));
            Zero = zeros(size(substructure.m.BB,1), size(phi_bar.N,2));
            alpha = [I, Zero; phi_bar.C, phi_bar.N];
        end

        function [p] = findSubstructureGeneralizedCoordinates(obj, substructure, alpha)
            % Compute the generalized coordinates of the substucture p from
            % its physical coodinates delta
            u = [substructure.u.B; substructure.u.I];
            p = alpha \ u;
        end

        function [m_bar, k_bar, P] = computeSubstructureGeneralizedMatrices(obj, substructure, phi_bar, w, alpha)
            % Compute mass and stiffness matrices, and force vector of the 
            % substructure in the generalized coordinates p
            m_bar = obj.computeSubstructureGeneralizedMass(substructure.m, alpha, phi_bar);
            k_bar = obj.computeSubstructureGeneralizedStiffness(substructure.k, alpha, phi_bar);
            P = obj.computeSubstructureGeneralizedForce(substructure.F, alpha);
            % Verify mass and stiffness using relation: kj_bar_NN = wj^2 * mj_bar_NN 
            k_bar_NN = zeros(length(w));
            for j = 1:length(w)
                k_bar_NN(j,j) = w(j)^2 * m_bar.NN(j,j);
            end
            assert(norm(k_bar.NN - k_bar_NN, 'fro') < obj.tolerance);
            % Final check to ensure correct sizes
            mbar = [m_bar.BB, m_bar.BN; m_bar.NB, m_bar.NN];
            kbar = [k_bar.BB, k_bar.BN; k_bar.NB, k_bar.NN];
            assert(all(size(mbar) == size(kbar)) && all([size(kbar, 1), 1] == size(P)));
        end

        function [m_bar] = computeSubstructureGeneralizedMass(obj, mr, alpha, phi_bar)
            % Compute mass matrix in the generalized coordinates p
            mBB = mr.BB + phi_bar.C' * mr.II * phi_bar.C;
            mBN = phi_bar.C' * mr.II * phi_bar.N;
            mNB = mBN';
            mNN = phi_bar.N' * mr.II * phi_bar.N;
            mbar = [mBB, mBN; mNB, mNN];
            m = [mr.BB, zeros(size(mr.BB,1), size(mr.II,2)); zeros(size(mr.II,1), size(mr.BB,2)), mr.II];
            m_bar = alpha' * m * alpha;  
            assert(norm(mbar - m_bar, 'fro') < obj.tolerance);
            m_bar = struct('BB', mBB, 'BN', mBN, 'NB', mNB, 'NN', mNN);       
        end

        function [k_bar] = computeSubstructureGeneralizedStiffness(obj, kr, alpha, phi_bar)
            % Compute stiffness matrix in the generalized coordinates p
            kBB = kr.BB + kr.BI * phi_bar.C;
            kNN = phi_bar.N' * kr.II * phi_bar.N;
            kBN = zeros(size(kBB, 1), size(kNN, 2));
            kNB = zeros(size(kNN, 1), size(kBB, 2));
            kbar = [kBB, kBN; kNB, kNN];
            k = [kr.BB, kr.BI; kr.IB, kr.II]; 
            k_bar = alpha' * k * alpha;
            assert(norm(kbar - k_bar, 'fro') < obj.tolerance);
            k_bar = struct('BB', kBB, 'BN', kBN, 'NB', kNB, 'NN', kNN);
        end

        function [P] = computeSubstructureGeneralizedForce(obj, f, alpha)
            % Compute force vector in the generalized coordinates
            F = [f.B; f.I];
            P = alpha' * F;
        end

        function [obj] = assembleStructureGlobalMatrices(obj)
            % Assemble mass and stiffness matrices and force vector of the structure
            % Create the connectivity matrix of the structures
            fprintf('%s Assembling generalized global matrices of the structure\n', repmat('-', 1, 8));
            i = 1;
            indices = zeros(numEntries(obj.substructures), 2);
            M = zeros(length(obj.structure.M));
            K = zeros(length(obj.structure.K));
            for r = 1:numEntries(obj.substructures)
                Sr = obj.substructures(r);
                m_bar = [Sr.m_bar.BB, Sr.m_bar.BN; Sr.m_bar.NB, Sr.m_bar.NN];
                k_bar = [Sr.k_bar.BB, Sr.k_bar.BN; Sr.k_bar.NB, Sr.k_bar.NN];
                N = size(m_bar, 1);
                M(i:i+N-1, i:i+N-1) = m_bar;
                K(i:i+N-1, i:i+N-1) = k_bar;
                indices(r, :) = [i, i + N - 1];
                i = i + N;
            end
            [obj.structure.M, obj.structure.K] = obj.enforceSubstructurBoundariesCompatibility(M, K, indices);
        end

        function [M, K] = enforceSubstructurBoundariesCompatibility(obj, M, K, indices)
            % Enforce geometric compatibility (continuity) between substructure boundaries
            fprintf('%s Enforcing continuity between substructures\n', repmat('-', 1, 10));
            for i = 1:numEntries(obj.substructures)
                for j = i+1:numEntries(obj.substructures)
                    if isequal(obj.substructures(i).nodes.Bd, obj.substructures(j).nodes.Bd)
                        % Enforce continuity for shared DOFs
                        uli = obj.substructures(i).dofs.l.Bd;
                        ulj = obj.substructures(j).dofs.l.Bg;
                        for k = 1:length(uli)
                            % Shared DoF in global matrices
                            ugi = indices(i, 1) + uli(k);
                            ugj = indices(j, 1) + ulj(k);

                            % Modify global mass and stiffness matrices to couple shared DOFs
                            M(ugi, :) = M(ugi, :) + M(ugj, :);
                            M(ugj, :) = M(ugi, :);
                            M(:, ugi) = M(:, ugi) + M(:, ugj);
                            M(:, ugj) = M(:, ugi);

                            K(ugi, :) = K(ugi, :) + K(ugj, :);
                            K(ugj, :) = K(ugi, :);
                            K(:, ugi) = K(:, ugi) + K(:, ugj);
                            K(:, ugj) = K(:, ugi);
                        end
                    end
                end
            end
        end

        function [obj] = computeStructureModalProperties(obj)
            % Compute the natural frequencies and modes shapes of the
            % structures using the reduced global matrices
            fprintf('%s Finding modes and mode shapes of the structure\n', repmat('-', 1, 6));
            % Finding and sort the eigenvalues and eigenvectors for free nodes
            [phi, Lambda] = eig(obj.structure.K, obj.structure.M);
            [w2, idx] = sort(diag(Lambda));
            obj.Phi = phi(:, idx);
            obj.fn = real(sqrt(w2) / 2*pi);
        end
    end
end
