%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary Function - Solve 2D thin plate modal analysis using ACM element
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef PlateModel
    % PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        % Geometry
        h
        H
        L
        Iz % cross section moment of inertia per length [m^3]
        
        % Material
        E
        nu
        rho
        
        % Mesh Definition
        eSize
        nElementX
        nElementZ
        nElements
        dx
        dz
        nNodesX
        nNodesZ
        nNodes
        ndof
        node_dof % Dof associated with each node
        node_coor % Node coordinates dataframe
        element_con  % Element conectivity dataframe
        index_table % DoF index dataframe
        
        % System matrices
        m_e % Element mass matrix
        k_e % Element stiffness
        M  % Global mass matrix
        K  % Global stiffness matrix
        
        % Boundary conditions
        fixed_nodes  % List of fixed nodes
        fixed_dof  % List of restrained DoFs
        
        % Results
        results
    end
    
    methods (Static)
        function obj = PlateModel(parameters) % Class Constructor
            % Geometry
            obj.h = parameters('h');
            obj.H = parameters('H');
            obj.L = parameters('L');
            obj.Iz = (obj.h ^ 3) / 12;
            % Material
            obj.E = parameters('E');
            obj.nu = parameters('nu');
            obj.rho = parameters('rho');
            % Model
            obj.eSize = parameters('eSize');
            
            % Mesh definitions
            obj.nElementX = uint16(obj.L / obj.eSize);
            obj.nElementZ = uint16(obj.H / obj.eSize);
            obj.nElements = obj.nElementX * obj.nElementZ;
            
            obj.dx = obj.L / double(obj.nElementX);
            obj.dz = obj.H / double(obj.nElementZ);
            
            obj.nNodesX = obj.nElementX + 1;
            obj.nNodesZ = obj.nElementZ + 1;
            obj.nNodes = obj.nNodesX * obj.nNodesZ;
            
            obj.ndof = obj.nNodes * 3; % w, theta_x and theta_y
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = get_element_matrices(obj)
            % Element mass matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            a = obj.dx / 2;
            b = obj.dz / 2;
            m11 = [3454    922*b     -922*a     1226     398*b     548*a
                922*b   320*b^2   -252*a*b   398*b    160*b^2   168*a*b
                -922*a  -252*a*b    320*a^2  -548*a   -168*a*b  -240*a^2
                1226    398*b     -548*a     3454     922*b     922*a
                398*b   160*b^2   -168*a*b   922*b    320*b^2   252*a*b
                548*a   168*a*b   -240*a^2   922*a    252*a*b   320*a^2];
            
            m21 =[394     232*b     -232*a     1226     548*b     398*a
                -232*b  -120*b^2    112*a*b  -548*b   -240*b^2  -168*a*b
                232*a   112*a*b   -120*a^2   398*a    168*a*b   160*a^2
                1226    548*b     -398*a     394      232*b     232*a
                -548*b  -240*b^2    168*a*b  -232*b   -120*b^2  -112*a*b
                -398*a  -168*a*b    160*a^2  -232*a   -112*a*b  -120*a^2];
            
            m22 = [3454   -922*b      922*a     1226    -398*b    -548*a
                -922*b   320*b^2   -252*a*b  -398*b    160*b^2   168*a*b
                922*a  -252*a*b    320*a^2   548*a   -168*a*b  -240*a^2
                1226   -398*b      548*a     3454    -922*b    -922*a
                -398*b   160*b^2   -168*a*b  -922*b    320*b^2   252*a*b
                -548*a   168*a*b   -240*a^2  -922*a    252*a*b   320*a^2];
            
            m = [m11 m21'
                m21 m22];
            
            obj.m_e =  (obj.rho * obj.h * a * b / 6300) * m;
            
            % Element stiffness Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            alf = obj.dx / obj.dz;
            bet = obj.dz / obj.dx;
            nu_ = obj.nu;
            I_1 = double(diag([-1  1  1]));
            I_2 = double(diag([ 1 -1  1]));
            I_3 = double(diag([ 1  1 -1]));
            
            k11 = [4*(bet^2+alf^2)+(2/5)*(7-2*nu_)      2*(2*alf^2+(1/5)*(1+4*nu_))*b          2*(-2*bet^2-(1/5)*(1+4*nu_))*a
                2*(2*alf^2+(1/5)*(1+4*nu_))*b           4*((4/3)*alf^2+(4/15)*(1-nu_))*b^2    -4*nu_*a*b
                2*(-2*bet^2-(1/5)*(1+4*nu_))*a         -4*nu_*a*b                              4*((4/3)*bet^2+(4/15)*(1-nu_))*a^2];
            
            k21 = [-2*(2*bet^2-alf^2)-(2/5)*(7-2*nu_)   2*(alf^2-(1/5)*(1+4*nu_))*b            2*(2*bet^2+(1/5)*(1-nu_))*a
                2*(alf^2-(1/5)*(1+4*nu_))*b             4*((2/3)*alf^2-(4/15)*(1-nu_))*b^2     0
                -2*(2*bet^2+(1/5)*(1-nu_))*a            0                                      4*((2/3)*bet^2-(1/15)*(1-nu_))*a^2];
            
            k31 = [-2*(bet^2+alf^2)+(2/5)*(7-2*nu_)     2*(-alf^2+(1/5)*(1-nu_))*b             2*(bet^2-(1/5)*(1-nu_))*a
                2*(alf^2-(1/5)*(1-nu_))*b               4*((1/3)*alf^2+(1/15)*(1-nu_))*b^2     0
                2*(-bet^2+(1/5)*(1-nu_))*a              0                                      4*((1/3)*bet^2+(1/15)*(1-nu_))*a^2];
            
            k41 = [2*(bet^2-2*alf^2)-(2/5)*(7-2*nu_)    2*(-2*alf^2-(1/5)*(1-nu_))*b           2*(-bet^2+(1/5)*(1+4*nu_))*a
                2*(2*alf^2+(1/5)*(1-nu_))*b             4*((2/3)*alf^2-(1/15)*(1-nu_))*b^2     0
                2*(-bet^2+(1/5)*(1+4*nu_))*a            0                                      4*((2/3)*bet^2-(4/15)*(1-nu_))*a^2];
            
            
            k22 = I_3'*k11*I_3;
            k32 = I_3'*k41*I_3;
            k42 = I_3'*k31*I_3;
            
            k33 = I_1'*k11*I_1;
            k43 = I_1'*k21*I_1;
            
            k44 = I_2'*k11*I_2;
            
            k = [k11  k21'  k31'  k41'
                k21  k22   k32'  k42'
                k31  k32   k33   k43'
                k41  k42   k43   k44];
            
            obj.k_e = ((obj.E*obj.h^3)/(48*(1 - obj.nu^2)*a*b)) * k;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = generate_mesh(obj)
            fprintf('...Generating mesh\n')
            % Generate the table of DoF IDs and corresponding nodes
            dof_per_node = 3;
            nodes = 1:obj.nNodes;
            dof = 1:(obj.nNodes*dof_per_node);
            dof = reshape(dof,dof_per_node,[]).';
            obj.node_dof = [nodes' dof];
            
            % Nodes coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.node_coor = zeros(obj.nNodes, 2);
            n = 1;
            for j = 0:obj.nNodesZ-1
                j = double(j);
                for i = 0:obj.nNodesX-1
                    i = double(i);
                    obj.node_coor(n, 1) = i * obj.dx;
                    obj.node_coor(n, 2) = j * obj.dz;
                    n = n + 1;
                end
            end
            
            % Conectivity matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Relates local nodes of elements to global nodes
            element_con_aux = zeros(obj.nElements, 4);
            
            n=1;
            m=1;
            for i=1:obj.nElementZ
                for j=1:obj.nElementX
                    element_con_aux(n,1) = m;
                    element_con_aux(n,2) = m + 1;
                    element_con_aux(n,3) = m + obj.nNodesX + 1 ;%nelx + 2;
                    element_con_aux(n,4) = m + obj.nNodesX;% + 1;
                    m=m+1;
                    n=n+1;
                end
                m = i*obj.nNodesX + 1;
            end
            
            obj.element_con = element_con_aux;
            clear element_con_aux n q r a
            % Index table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Relates the gloobal DoF associated with each element
            index_table_aux = zeros(obj.nElements, 12);
            for n = 1:obj.nElements
                index_table_aux(n, 1) = obj.node_dof(obj.element_con(n, 1), 2); % w
                index_table_aux(n, 2) = obj.node_dof(obj.element_con(n, 1), 3); % thetaX
                index_table_aux(n, 3) = obj.node_dof(obj.element_con(n, 1), 4); % thetaZ
                
                index_table_aux(n, 4) = obj.node_dof(obj.element_con(n, 2), 2);
                index_table_aux(n, 5) = obj.node_dof(obj.element_con(n, 2), 3);
                index_table_aux(n, 6) = obj.node_dof(obj.element_con(n, 2), 4);
                
                index_table_aux(n, 7) = obj.node_dof(obj.element_con(n, 3), 2);
                index_table_aux(n, 8) = obj.node_dof(obj.element_con(n, 3), 3);
                index_table_aux(n, 9) = obj.node_dof(obj.element_con(n, 3), 4);
                
                index_table_aux(n, 10) = obj.node_dof(obj.element_con(n, 4), 2);
                index_table_aux(n, 11) = obj.node_dof(obj.element_con(n, 4), 3);
                index_table_aux(n, 12) = obj.node_dof(obj.element_con(n, 4), 4);
            end
            obj.index_table = index_table_aux;
            clear index_table_aux n
            
            
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = get_global_matrices(obj)
            fprintf('...Getting global matrices\n');
            K_aux = zeros(obj.ndof, obj.ndof);
            M_aux = zeros(obj.ndof, obj.ndof);
            
            for n = 1:obj.nElements
                idx = obj.index_table(n,:); % Grid of indexes
                K_aux(idx, idx) = K_aux(idx, idx) + obj.k_e;
                M_aux(idx, idx) = M_aux(idx, idx) + obj.m_e;
            end
            
            obj.K = sparse(K_aux);
            obj.M = sparse(M_aux);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = apply_bc(obj, regions)
            % Remove fixed DoFs from global matrices
            %
            % Args:
            %     regions ((1,4) array: coordinates of Two Points Box
            %     as in {n: [[x1,z1],[x2,z2]])}. All DoF of nodes
            %     inside the box will be removed from the global system of equations.
            fprintf('...Applying BC to global matrices\n')
            tol = [obj.dx/2 obj.dz/2];
            obj.fixed_nodes = [];
            for k = keys(regions)
                r = regions(k{1});
                idx = find(obj.node_coor(:,1) >= r(1) - tol(1) & ...
                    obj.node_coor(:,1) <= r(3) + tol(1) & ...
                    obj.node_coor(:,2) >= r(2) - tol(2) & ...
                    obj.node_coor(:,2) <= r(4) + tol(2));
                obj.fixed_nodes = [obj.fixed_nodes idx.'];
            end
            
            obj.fixed_dof = [];
            for node = obj.fixed_nodes
                obj.fixed_dof = [obj.fixed_dof obj.node_dof(node, :)]; % concatenate
            end
            obj.fixed_dof = reshape(obj.fixed_dof,1,[]);
            obj.fixed_dof = sort(obj.fixed_dof);
            
            % Removing fixed DoF from [Q] and [H]
            mask = true(obj.ndof,1);
            mask(obj.fixed_dof) = false;
            obj.M = obj.M(mask,:);
            obj.M = obj.M(:,mask);
            obj.K = obj.K(mask,:);
            obj.K = obj.K(:,mask);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = solve_eigenvalue_problem(obj, N_modes)
            fprintf('...Solving eigenvalue problem\n');
            [Vc, Wn2] = eigs(obj.K, obj.M, N_modes, 'sm');
            fn = diag(Wn2).^(0.5)/(2*pi); % in Hertz
            % Ordering eigenvalues and eigenvectors
            [fn, idx] = sort(fn);
            Vc = Vc(:, idx);
            
            % Normalizing eigenvectors matrix by the mass matrix,
            % such that Vc' * M * Vc = I
            m_r = diag(Vc'*obj.M*Vc);
            m_r = 1./(m_r.^0.5);
            for i = 1:size(Vc,2)
                Vc(:,i) = Vc(:,i).*m_r(i);
            end
            
            
            % Adding fixed Dof back into the displacement vector
            [~, nModos] = size(Vc);
            fixed = zeros(1,nModos);
            for row = obj.fixed_dof                
                Vc = [Vc(1:row-1,:) ; fixed ; Vc(row:end,:)];
            end
            
            results_aux.fn = real(fn);
            results_aux.Vc = real(Vc(1:3:end,:)); % Displacement shape (dof 1,4,7,10,13,...)
            
            obj.results = results_aux;
        end
                
    end
end