%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary Function - Solve 3D acoustic cavity modal analysis using
% regular hexa elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef AcousticModel
    % PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        W % volume width [m] (X axis)
        D% volume depth [m] (Y axis)
        H % volume height [m] (Z axis)
        % Fluid
        c % sound speed [m/s]
        rho % density [kg/m**3]
        % Model
        eSize % element size [m]
        % Mesh Definition
        nElementX
        nElementY
        nElementZ
        nElements
        dx
        dy
        dz
        nNodesX
        nNodesY
        nNodesZ
        nNodes
        ndof % Total number of DoF
        node_dof  % ID of Dof associated with each node
        node_coor  % Node coordinates dataframe
        element_con  % Element conectivity dataframe
        index_table  % DoF index dataframe
        % Element matrices
        H_e % Acoustic stiffness matrix
        Q_e % Acoustic inertia matrix
        % Global matrices
        Qg % Global acoustic inertia matrix (sparse)
        Hg % Global acoustic stiffness matrix (sparse)
        % Boundary conditions
        fixed_nodes % List of fixed nodes
        fixed_dof % List of restrained DoFs
        % Results
        results % Struct with eigenvalues and eigenvectors
    end
    
    % METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = AcousticModel(parameters) % Class Constructor
            % Geometry
            obj.W = parameters('W');  % volume width [m] (X axis)
            obj.D = parameters('D');  % volume depth [m] (Y axis)
            obj.H = parameters('H');  % volume height [m] (Z axis)
            
            % Fluid
            obj.c = parameters('c');  % sound speed [m/s]
            obj.rho = parameters('rho');  % density [kg/m**3]
            
            % Model
            obj.eSize = parameters('eSize');  % element size [m]
            
            % Mesh Definition
            obj.nElementX = uint16(obj.W / obj.eSize);
            obj.nElementY = uint16(obj.D / obj.eSize);
            obj.nElementZ = uint16(obj.H / obj.eSize);
            
            obj.nElements = obj.nElementX * obj.nElementY * obj.nElementZ;
            
            obj.dx = obj.W / double(obj.nElementX);
            obj.dy = obj.D / double(obj.nElementY);
            obj.dz = obj.H / double(obj.nElementZ);
            
            obj.nNodesX = obj.nElementX + 1;
            obj.nNodesY = obj.nElementY + 1;
            obj.nNodesZ = obj.nElementZ + 1;
            
            obj.nNodes = obj.nNodesX * obj.nNodesY * obj.nNodesZ;
            obj.ndof = obj.nNodes;  % Total number of DoF
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = get_element_matrices(obj)
            Be = [1 -1 -1 -1 1 1 1 -1;
                1 1 -1 -1 -1 -1 1 1;
                1 1 1 -1 1 -1 -1 -1;
                1 -1 1 -1 -1 1 -1 1;
                1 -1 -1 1 1 -1 -1 1;
                1 1 -1 1 -1 1 -1 1;
                1 1 1 1 1 1 1 1;
                1 -1 1 1 -1 -1 1 -1;];
            a1 = obj.dx / 2;
            a2 = obj.dy / 2;
            a3 = obj.dz / 2;
            
            % Eq. 8.89c,d,e (Fahy), but q is a column vector
            syms xi1 xi2 xi3
            q = [1; xi1; xi2; xi3; xi1*xi2; xi1*xi3; xi2*xi3; xi1*xi2*xi3];
            q1 = [0; 1; 0; 0; xi2; xi3; 0; xi2 * xi3];
            q2 = [0; 0; 1; 0; xi1; 0; xi3; xi1 * xi3];
            q3 = [0; 0; 0; 1; 0; xi1; xi2; xi1 * xi2];
            
            % Matrizes elementares
            obj.H_e = double(inv(Be')*int(int(int((((1/a1^2)*q1*q1') + ((1/a2^2)*q2*q2') + ((1/a3^2)*q3*q3'))*(a1*a2*a3),xi1,-1,1),xi2,-1,1),xi3,-1,1)*inv(Be));
            obj.Q_e = double((a1*a2*a3/obj.c^2)*inv(Be')*int(int(int(q*q',xi1,-1,1),xi2,-1,1),xi3,-1,1)*inv(Be));
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = generate_mesh(obj)
            fprintf('...Generating mesh\n')
            % Generate the table of DoF IDs and corresponding nodes
            dof_per_node = 1;
            nodes = 1:obj.nNodes;
            dof = 1:(obj.nNodes*dof_per_node);
            dof = reshape(dof,[],dof_per_node);
            obj.node_dof = [nodes' dof];
            
            % Nodes coordinates
            obj.node_coor = zeros(obj.nNodes, 3);
            n = 1;
            for k = 0:obj.nNodesZ-1
                k = double(k);
                for j = 0:obj.nNodesY-1
                    j = double(j);
                    for i = 0:obj.nNodesX-1
                        i = double(i);
                        obj.node_coor(n, 1) = i * obj.dx;
                        obj.node_coor(n, 2) = j * obj.dy;
                        obj.node_coor(n, 3) = k * obj.dz;
                        n = n + 1;
                    end
                end
            end
            
            
            % Conectivity matrix
            % Relates local nodes of elements to global nodes
            element_con_aux = zeros(obj.nElements, 8);
            n = 1;
            for k = 1:obj.nElementZ
                % k-th slice in the Z plane
                for j = 1:obj.nElementY
                    % j-th slice in the Y plane
                    m = (k-1) * obj.nNodesX * obj.nNodesY + (j-1) * obj.nNodesX + 1;
                    for i = 1:obj.nElementX
                        element_con_aux(n, 1) = m;
                        element_con_aux(n, 2) = m + 1;
                        element_con_aux(n, 3) = m + obj.nNodesX + 1;
                        element_con_aux(n, 4) = m + obj.nNodesX;
                        element_con_aux(n, 5) = m + obj.nNodesX * obj.nNodesY;
                        element_con_aux(n, 6) = m + obj.nNodesX * obj.nNodesY + 1;
                        element_con_aux(n, 7) = m + obj.nNodesX * obj.nNodesY + obj.nNodesX + 1;
                        element_con_aux(n, 8) = m + obj.nNodesX * obj.nNodesY + obj.nNodesX;
                        n = n + 1;
                        m = m + 1;
                    end
                end
            end
            obj.element_con = element_con_aux;
            clear element_con_aux i j k n m
            
            % Index table
            % Relates the gloobal DoF associated with each element
            index_table_aux = zeros(obj.nElements, 8);
            for n = 1:obj.nElements
                index_table_aux(n, 1) = obj.node_dof(obj.element_con(n, 1), 2);
                index_table_aux(n, 2) = obj.node_dof(obj.element_con(n, 2), 2);
                index_table_aux(n, 3) = obj.node_dof(obj.element_con(n, 3), 2);
                index_table_aux(n, 4) = obj.node_dof(obj.element_con(n, 4), 2);
                index_table_aux(n, 5) = obj.node_dof(obj.element_con(n, 5), 2);
                index_table_aux(n, 6) = obj.node_dof(obj.element_con(n, 6), 2);
                index_table_aux(n, 7) = obj.node_dof(obj.element_con(n, 7), 2);
                index_table_aux(n, 8) = obj.node_dof(obj.element_con(n, 8), 2);
            end
            obj.index_table = index_table_aux;
            clear index_table_aux n
            
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = get_global_matrices(obj)
            fprintf('...Getting global matrices\n');
            Q_aux = zeros(obj.ndof, obj.ndof);
            H_aux = zeros(obj.ndof, obj.ndof);
            
            for n = 1:obj.nElements
                idx = obj.index_table(n,:); % Grid of indexes
                Q_aux(idx, idx) = Q_aux(idx, idx) + obj.Q_e;
                H_aux(idx, idx) = H_aux(idx, idx) + obj.H_e;
            end
            
            obj.Qg = sparse(Q_aux);
            obj.Hg = sparse(H_aux);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = apply_bc(obj, regions)
            % Remove fixed DoFs from global matrices
            %
            % Args:
            %     regions ((1,4) array: coordinates of Two Points Box
            %     as in {n: [[x1,y1,z1],[x2,y2,z2]])}. All DoF of nodes
            %     inside the box will be removed from the global system of equations.
            fprintf('...Applying BC to global matrices\n')
            tol = [obj.dx/2 obj.dy/2 obj.dz/2];
            obj.fixed_nodes = [];
            for k = keys(regions)
                r = regions(k{1});
                idx = find(obj.node_coor(:,1) >= r(1,1) - tol(1) & ...
                           obj.node_coor(:,1) <= r(2,1) + tol(1) & ...
                           obj.node_coor(:,2) >= r(1,2) - tol(2) & ...
                           obj.node_coor(:,2) <= r(2,2) + tol(2) & ...
                           obj.node_coor(:,3) >= r(1,3) - tol(3) & ...
                           obj.node_coor(:,3) <= r(2,3) + tol(3));
                obj.fixed_nodes = [obj.fixed_nodes idx];
            end            
            obj.fixed_nodes = unique(obj.fixed_nodes);
            
            obj.fixed_dof = [];
            for node = obj.fixed_nodes
                obj.fixed_dof = [obj.fixed_dof obj.node_dof(node, :)]; % concatenate
            end
            obj.fixed_dof = reshape(obj.fixed_dof,1,[]);
            obj.fixed_dof = unique(obj.fixed_dof);
            obj.fixed_dof = sort(obj.fixed_dof);
            
            % Removing fixed DoF from [Q] and [H] using a boolean mask
            mask = true(obj.ndof, 1);
            mask(obj.fixed_dof) = false;
            obj.Qg = obj.Qg(mask,:);
            obj.Qg = obj.Qg(:,mask);
            obj.Hg = obj.Hg(mask,:);
            obj.Hg = obj.Hg(:,mask);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = solve_eigenvalue_problem(obj, N_modes)
            fprintf('...Solving eigenvalue problem\n');
            [Vc, Wn2] = eigs(obj.Hg, obj.Qg, N_modes, 'sm');
            fn = diag(Wn2).^(0.5)/(2*pi); % in Hertz
            % Ordering eigenvalues and eigenvectors
            [fn, idx] = sort(fn);
            Vc = Vc(:, idx);
            
            % Normalizing eigenvectors matrix by the mass matrix,
            % such that Vc' * M * Vc = I
            m_r = diag(Vc'*obj.Qg*Vc);
            m_r = 1./(m_r.^0.5);
            for i = 1:size(Vc,2)
                Vc(:,i) = Vc(:,i).*m_r(i);
            end
            
            % Adding fixed Dof back into the displacement vector
            if ~isempty(obj.fixed_dof)
                [~, nModos] = size(Vc);
                fixed = zeros(1,nModos);
                for row = obj.fixed_dof
                    Vc = [Vc(1:row-1,:) ; fixed ; Vc(row:end,:)];
                end
            end
            
            results_aux.fn = real(fn);
            results_aux.Phi = real(Vc);
            obj.results = results_aux;
        end
        
    end
end