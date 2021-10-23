close all; clc; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 1 - Calcule as 5 primeiras frequ^encias naturais e formas modais da
% placa e da cavidade quando desacopladas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Part 1 - Plate mode shape
fprintf('Solving door model\n')
door_param = containers.Map;
door_param('h') = 1.6e-3;  % Thickness [m]
door_param('H') = 500e-3;  % Door height [m]
door_param('L') = 200e-3;  % Door length [m]
door_param('E') = 186e9;   % Modulus of elasticity [Pa]
door_param('nu') = 0.3;    % Poison coef
door_param('rho') = 7870;  % Steel density [kg/m**3]
door_param('eSize') = 0.02; % Element size [m]

door = PlateModel(door_param);
door = door.get_element_matrices(door);
door = door.generate_mesh(door);
door = door.get_global_matrices(door);

% Fixed regions (Two-points recatangles [x1 z1, x2 z2])
regions=containers.Map;
regions('1')=[0.000 0.100, 0.020 0.120]; % Lower hinge
regions('2')=[0.000 0.380, 0.020 0.400]; % Upper hinge
regions('3')=[0.180 0.260, 0.220 0.280]; % Knob

door = door.apply_bc(door, regions);
door = door.solve_eigenvalue_problem(door, 50);

%% 1.1 - Plot Door mode shapes
x = 0 : door.dx : door.L;
z = 0 : door.dz : door.H;
[X, Z] = meshgrid(x, z);
for mode=1:5
    figure()
    mode_shape = door.results.Vc(:,mode);
    mode_shape = reshape(mode_shape, door.nNodesX, door.nNodesZ);
    surf(X, Z, mode_shape.');
    grid on
    box on;
end

%% Part 2 - Solve acoutic cavity model
fprintf('Solving acoustic cavity model\n')

cvt_param = containers.Map;
cvt_param('W') = 240e-3; % volume width [m] (X axis
cvt_param('D') = 240e-3; % volume depth [m] (Y axis)
cvt_param('H') = 540e-3; % volume height [m] (Z axis)
cvt_param('c') = 348; % sound speed [m/s]
cvt_param('rho') = 1.7; % density [kg/m**3]
cvt_param('eSize') = 0.02; % element size [m]

cvt = AcousticModel(cvt_param);
cvt = cvt.get_element_matrices(cvt);
cvt = cvt.generate_mesh(cvt);
cvt = cvt.get_global_matrices(cvt);
cvt = cvt.solve_eigenvalue_problem(cvt, 5);


%% 2.1 - Plot decoupled acoustic cavity mode shapes
Vc = cvt.results.Phi;
x = 0 : cvt.dx : cvt.W;
y = 0 : cvt.dy : cvt.D;
z = 0 : cvt.dz : cvt.H;
[X, Y, Z] = meshgrid(x, y, z);

for mode = 1:5    
    figure()
    set(gcf, 'Position', get(0, 'Screensize'));
    mode_shape = reshape(Vc(:,mode), [cvt.nNodesX, cvt.nNodesY, cvt.nNodesZ]);
    mode_shape = permute(mode_shape,[2 1 3]);
    slice(X,Y,Z, mode_shape, 0, 0, 0.540);
    box on
    ax = gca;
    ax.BoxStyle = 'full';
    shading interp
    set(gca, 'FontSize', 20)
    xlabel('Width [m]')
    ylabel('Length [m]')
    zlabel('Height [m]')
    title(sprintf('Decoupled acoustic cavity - Mode %.f (%0.1f [Hz])',mode,cvt.results.fn(mode)));
    colormap jet;colorbar
    set(gcf,'color','w');
    axis equal
end

%% Part 3 - Solve acoutic cavity with back opening
fprintf('Solve acoutic cavity with back opening\n')

cvt2 = AcousticModel(cvt_param);
cvt2 = cvt2.get_element_matrices(cvt2);
cvt2 = cvt2.generate_mesh(cvt2);
cvt2 = cvt2.get_global_matrices(cvt2);

% Constant pressure regions (Two-points box [x1 y1 z1, x2 y2 z2])
regions = containers.Map;
regions('1')=[0.180 0.240 0.040; 0.220 0.240 0.100]; % back opening

cvt2 = cvt2.apply_bc(cvt2, regions);
cvt2 = cvt2.solve_eigenvalue_problem(cvt2, 5);

%% 3.1 - Plot
Vc = cvt2.results.Phi;
x = 0 : cvt.dx : cvt.W;
y = 0 : cvt.dy : cvt.D;
z = 0 : cvt.dz : cvt.H;
[X, Y, Z] = meshgrid(x, y, z);

for mode = 1:5    
    figure()
    set(gcf, 'Position', get(0, 'Screensize'));
    mode_shape = reshape(Vc(:,mode), [cvt.nNodesX, cvt.nNodesY, cvt.nNodesZ]);
    mode_shape = permute(mode_shape,[2 1 3]);
    xslice = [0.000, 0.240];   
    yslice = [0.000, 0.240];
    zslice = [0.000, 0.540];
    slice(X,Y,Z, mode_shape, xslice, yslice, zslice);
    box on;
    grid on;
    ax = gca;
    ax.BoxStyle = 'full';
    shading interp
    set(gca, 'FontSize', 20)
    xlabel('Width [m]')
    ylabel('Length [m]')
    zlabel('Height [m]')
    title(sprintf('Decoupled acoustic cavity - Mode %.f (%0.1f [Hz])',mode,cvt.results.fn(mode)));
    colormap jet;colorbar
    set(gcf,'color','w');
    axis equal
end

%% Part 4 - Strong-cloupling
a = door.dx;
b = door.dz;
A_e = [1  -1  -1   1   1    1  -1   -1   -1  -1    1    1 ;
       0   0  1/b  0 -1/b -2/b  0   1/b  2/b 3/b -1/b -3/b;
       0 -1/a  0  2/a 1/a   0 -3/a -2/a -1/a  0   3/a  1/a;
       1   1  -1   1  -1    1   1   -1    1  -1   -1   -1 ;
       0   0  1/b  0  1/b -2/b  0   1/b -2/b 3/b  1/b  3/b;
       0 -1/a  0 -2/a 1/a   0 -3/a  2/a -1/a  0   3/a  1/a;
       1   1   1   1   1    1   1    1    1   1    1    1 ;
       0   0  1/b  0  1/b  2/b  0   1/b  2/b 3/b  1/b  3/b;
       0 -1/a  0 -2/a -1/a  0 -3/a -2/a -1/a  0  -3/a -1/a;
       1  -1   1   1   -1   1  -1    1   -1   1   -1   -1 ;
       0   0  1/b  0  -1/b 2/b  0   1/b -2/b 3/b -1/b -3/b;
       0 -1/a  0  2/a -1/a  0 -3/a  2/a -1/a  0  -3/a -1/a];
   
B_e = [1 -1 -1 -1  1  1  1 -1;
       1  1 -1 -1 -1 -1  1  1;
       1  1  1 -1  1 -1 -1 -1;
       1 -1  1 -1 -1  1 -1  1;
       1 -1 -1  1  1 -1 -1  1;
       1  1 -1  1 -1  1 -1  1;
       1  1  1  1  1  1  1  1;
       1 -1  1  1 -1 -1  1 -1;];

% ACM element form functions
syms xi eta
p = [1; xi; eta; xi^2; xi*eta; eta^2; xi^3; xi^2*eta; xi*eta^2; eta^3; xi^3*eta; xi*eta^3];
% Hexa element form functions
syms xi1 xi2 xi3 
q = [1; xi1; xi2; xi3; xi1*xi2; xi1*xi3; xi2*xi3; xi1*xi2*xi3];
% Replacing the hexa (acoustic) coordinate system by the ACM
xi1 = xi;  
xi2 = eta;
xi3 = -1;
q = subs(q);
% Computing Element Coumpling Matrix (12x8)
S_e = double(inv(A_e')*int(int(p*q'*a*b,xi,-1,1),eta,-1,1)*inv(B_e));

% Assembling Global Coupling Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = zeros(door.ndof, cvt.ndof); % ndof: Number of DoF
% The door and the cavity meshes are not coincident, so we must perform a
% translation. Using cavity mesh coordinate system (csys) as the default.
door_cvt_nodes = zeros(door.nNodes,1);
tol = 0.005;
x_offset = 0.02;
y_offset = 0.02;
for node = 1:door.nNodes
    condition_1 = abs(cvt.node_coor(:,1) - (door.node_coor(node,1) + x_offset)) < tol;
    condition_2 = abs(cvt.node_coor(:,2)) < tol;
    condition_3 = abs(cvt.node_coor(:,3) - (door.node_coor(node,2) + y_offset)) < tol;
    idx = condition_1 & condition_2 & condition_3;
    if sum(idx)~=1
        fprintf('Meshes are not concident!!')
    else
        door_cvt_nodes(node) = find(idx==1); % Door node in cavity csys
    end
end
door_cvt_elements = zeros(door.nElements, 1);
for e = 1:door.nElements    
    n1 = cvt.element_con(:,1) == door_cvt_nodes(door.element_con(e,1));
    n2 = cvt.element_con(:,2) == door_cvt_nodes(door.element_con(e,2));
    n3 = cvt.element_con(:,6) == door_cvt_nodes(door.element_con(e,3));
    n4 = cvt.element_con(:,5) == door_cvt_nodes(door.element_con(e,4));
    idx = n1 & n2 & n3 & n4;
    if sum(idx)~=1
        fprintf('Meshes are not concident!!')
    else
        door_cvt_elements(e) = find(idx==1); % Door element in cavity csys
    end
end

% Coupling the I-th dof from door to the J-th dof from cavity
for e = 1:length(door.element_con) % For each door element
    z = 1;
    for j = 1:4 % For each node of the ACM element
        J = door.element_con(e, j); % Global node (in door csys)  
        for i = 1:3 
            % LM_ES(z) is the i-th door DoF in the door csys
            door_dof(z) = door.node_dof(J, i);
            z = z + 1;
        end          
    end
    
    r = 1; 
    e = door_cvt_elements(e); % Convert the node from door csys to cavity csys
    for j = 1:8 % For each node of the Hex element
        J = cvt.element_con(e, j);  % Global node (in cavity csys) 
        for i = 1
            % LM_AC(z) is the i-th cavity DoF in the cavity csys
            cvt_dof(r) = cvt.node_dof(J, i);
            r = r + 1;
        end
    end
    
    % Sum the node contribution to the global coupling matrix
    for i = 1:12
        I = door_dof(i); % i-th door DoF in door csys
        if I~=0
            for j = 1:8
                J = cvt_dof(j); % j-th cavity DoF in cavity csys
                if J~=0
                    S(I,J) = S(I,J) + S_e(i,j);
                end
            end
        end
    end
end

% Removing fixed DoF
mask = true(door.ndof,1); % Boolean mask for the rows
mask(door.fixed_dof) = false;
S = S(mask,:);
% Assembling the coupled matrices
M = door.M;
K = door.K;
Q = cvt.Qg;
H = cvt.Hg;
zeros1 = zeros(length(M),length(Q));
zeros2 = zeros(length(H),length(K));

M_coupled = [   M   zeros1 ;
               -S'    Q   ];
M_coupled = sparse(M_coupled);

K_coupled = [   K      S   ;
              zeros2   H  ];
K_coupled = sparse(K_coupled);

%% 4.1 - Solving the eigenvalue problem
[Vc, Wn2] = eigs(K_coupled, M_coupled, 10, 'sm');
fn = diag(Wn2).^(0.5)/(2*pi); % in Hertz
% Ordering eigenvalues and eigenvectors
[fn, idx] = sort(fn);
Vc = Vc(:, idx);
% Ignore complex part
fn = real(fn);
Vc = real(Vc);

% Normalizing eigenvectors matrix by the mass matrix,
% such that Vc' * M_coupled * Vc = I
m_r = diag(Vc'*M_coupled*Vc);
m_r = 1./(m_r.^0.5);
for i = 1:size(Vc,2)
    Vc(:,i) = Vc(:,i).*m_r(i);
end

%% 4.2 - Splitting the DoF vector
Vc_door = Vc(1:length(door.M), :);
% Adding fixed Dof back into the displacement vector
[~, nModos] = size(Vc_door);
fixed = zeros(1, nModos);
for row = door.fixed_dof
    Vc_door = [Vc_door(1:row-1,:) ; fixed ; Vc_door(row:end,:)];
end
% Taking only the vertical displacements DoF from door
Vc_door = real(Vc_door(1:3:end,:)); % Displacement shape (dof 1,4,7,10,13,...)

% Taking the cavity results
Vc_cvt = Vc(length(door.M)+1:end, :);

%% 4.3 - Plot Door mode shapes
x = 0 : door.dx : door.L;
z = 0 : door.dz : door.H;
[X, Z] = meshgrid(x, z);
for mode = 1:5
    figure()
    mode_shape = Vc_door(:,mode);
    mode_shape = reshape(mode_shape, door.nNodesX, door.nNodesZ);
    surf(X, Z, mode_shape.');
    grid on;
    box on;
end

%% 4.4 - Plot Cavity mode shapes
x = 0 : cvt.dx : cvt.W;
y = 0 : cvt.dy : cvt.D;
z = 0 : cvt.dz : cvt.H;
[X, Y, Z] = meshgrid(x, y, z);

for mode = 1:5    
    figure()
    set(gcf, 'Position', get(0, 'Screensize'));
    mode_shape = reshape(Vc_cvt(:,mode), [cvt.nNodesX, cvt.nNodesY, cvt.nNodesZ]);
    mode_shape = permute(mode_shape,[2 1 3]);
    slice(X,Y,Z, mode_shape, 0, 0, 0.540);
    
    box on
    ax = gca;
    ax.BoxStyle = 'full';
    shading interp
    set(gca, 'FontSize', 20)
    xlabel('Width [m]')
    ylabel('Length [m]')
    zlabel('Height [m]')
    title(sprintf('Decoupled acoustic cavity - Mode %.f (%0.1f [Hz])',mode,cvt.results.fn(mode)));
    colormap jet;colorbar
    set(gcf,'color','w');
    axis equal
end

%% 4.5 - Compute FRF
% Find nodes and DoF of interest
A_coor = [0.160 0.420]; % Force at the door. Door csys [x z]
B_coor = [0.100 0.160 0.340]; % Response inside cavity. Cavity csys [x y z]

x = door.node_coor(:,1) == A_coor(1);
z = door.node_coor(:,2) == A_coor(2);
idx = x & z;
if sum(idx)==1
    A_node = find(idx==1);
else
    fprintf('There is no node on the selected location')
end
A_dof = door.node_dof(A_node,2); % In door csys (also in coupled-system csys)

x = cvt.node_coor(:,1) == B_coor(1);
y = cvt.node_coor(:,2) == B_coor(2);
z = cvt.node_coor(:,3) == B_coor(3);
idx = x & y & z;
if sum(idx)==1
    B_node = find(idx==1);
else
    fprintf('There is no node on the selected location')
end
B_dof = cvt.node_dof(B_node,2); % In cavity csys
B_dof_coupled_csys = B_dof + length(M); % In coupled-system csys

% % FRF (m/N)
% % coupled
% % only one response point
% eta=0.03;
% freq = 10:0.5:500;
% H_jk = zeros(door.ndof/3, length(freq));
% k = find(door.node_dof(:,2)==A_dof); % DoF index without rotations
% for j = 1:door.ndof/3 % For each response dof
%     for p = 1:length(freq) % At a given frequency
%         for n = 1:length(fn) % Sum the contribution of every mode
%             H_jk(j, p) = H_jk(j, p) + Vc_door(j, n)*Vc_door(k, n)/...
%                                     (fn(n)^2 - freq(p)^2 + 1j*eta*fn(n)^2);
%         end
%     end
% end

% FRF (m/N) Decoupled
eta=0.03;
freq = 10:0.5:500;
H_jk = zeros(door.ndof/3, length(freq));
k = find(door.node_dof(:,2)==A_dof); % DoF index without rotations
for j = 1:size(door.results.Vc, 1) % For each response dof
    for p = 1:length(freq) % At a given frequency
        for n = 1:length(door.results.fn) % Sum the contribution of every mode
            H_jk(j, p) = H_jk(j, p) + door.results.Vc(j, n)*door.results.Vc(k, n)/...
                                    (door.results.fn(n)^2 - freq(p)^2 + 1j*eta*door.results.fn(n)^2);
        end
    end
end

% FRF (Pa/(m^3/m^3)) Decoupled
% frf between point B and every other DoF of the cavity
eta=0.03;
freq = 10:0.5:500;
G_jk = zeros(cvt.ndof, length(freq));
k = find(cvt.node_dof(:,2) == B_dof); % DoF index without rotations
for j = 1:size(cvt.results.Phi, 1) % For each response dof
    for p = 1:length(freq) % At a given frequency
        for n = 1:length(cvt.results.fn) % Sum the contribution of every mode
            G_jk(j, p) = G_jk(j, p) + cvt.results.Phi(j, n)*cvt.results.Phi(k, n)/...
                                    (cvt.results.fn(n)^2 - freq(p)^2 + 1j*eta*cvt.results.fn(n)^2);
        end
    end
end

G_jk_door = G_jk(door_cvt_nodes, :); % frf between point B and the DoF at the door (Pa/(m^3/m^3))
vol_element = cvt.dx*cvt.dy*cvt.dz;
srf_element = cvt.dx*cvt.dz;
P_jk_door = G_jk_door*vol_element/srf_element; % Pa/m (Pa per unit displacement of element face.

I = zeros(size(H_jk));
for j = 1:size(P_jk_door,1)
    I(j,:) = H_jk(j,:) .* P_jk_door(j,:);
end

semilogy(sum(real(I), 1))

% A_jk = zeros(length(freq),1);
% for p = 1:length(freq)
%     for n = 1:length(fn)
%         A_jk(p) = A_jk(p) + Vc(A_dof, n)*Vc(B_dof, n)/(fn(n)^2 - freq(p)^2 + 1j*eta*fn(n)^2);
%     end
% end
% 
% semilogy(abs(A_jk))
%% 4.6 - Plot ansys and code
Ansys_pointA = xlsread('Ansys_pointA.xlsx'); % Meter/Newton
Ansys_pointA = [Ansys_pointA(:,2) Ansys_pointA(:,3).*exp(1j*deg2rad(Ansys_pointA(:,4)))];
Ansys_pointB = xlsread('Ansys_pointB.xlsx'); % Pascal(B)/Newton(A)
Ansys_pointB = [Ansys_pointB(:,2) Ansys_pointB(:,3).*exp(1j*deg2rad(Ansys_pointB(:,4)))];

figure(); hold on;
f = real(Ansys_pointB(:,1));
semilogy(f, abs(Ansys_pointB(:,2)))

semilogy(freq, abs(sum(real(I), 1))*1e-5)

%% Plot [NOT WORKING}

% for mode = 2:2
%     figure();    hold on;
%     % Plot Cavity mode shapes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     ax1 = axes;
%     x = 0 : cvt.dx : cvt.W;
%     y = 0 : cvt.dy : cvt.D;
%     z = 0 : cvt.dz : cvt.H;
%     [X, Y, Z] = meshgrid(x, y, z);    
%     mode_shape = reshape(Vc_cvt(:,mode), [cvt.nNodesX, cvt.nNodesY, cvt.nNodesZ]);
%     mode_shape = permute(mode_shape,[2 1 3]);
%     slice(X,Y,Z, mode_shape, 0, 0, 0.540);
%     
%     box on
%     ax = gca;
%     ax.BoxStyle = 'full';
%     shading interp
%     set(gca, 'FontSize', 20)
%     xlabel('Width [m]')
%     ylabel('Length [m]')
%     zlabel('Height [m]')
%     title(sprintf('Coupled acoustic cavity - Mode %.f (%0.1f [Hz])',mode,cvt.results.fn(mode)));
% %     colormap jet;
% %     colorbar
%     set(gcf,'color','w');
%     axis equal
%     view(ax1,[-69.5 18.123046875])
%     
%      
%     % Plot Door Mode Shapes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     ax2 = axes;
%     set(gcf, 'Position', get(0, 'Screensize'));
%     x = 0 : door.dx : door.L;
%     z = 0 : door.dz : door.H;
%     x = x + door.dx; z = z + door.dz;
%     [X, Z] = meshgrid(x, z);
%     mode_shape = Vc_door(:,mode);
%     mode_shape = 0.04*mode_shape/max(abs(mode_shape));
%     mode_shape = reshape(mode_shape, door.nNodesX, door.nNodesZ);
%     surf(X, mode_shape.' - 0.1, Z);
%     grid on;
%     box on;
%     ax = gca;
%     ax.BoxStyle = 'full';
%     axis off;    
%     
%     
%     % Link two axes together
%     hLink = linkprop([ax1,ax2],{'CameraUpVector','CameraPosition','CameraTarget'});
%     % Hide the top axes
%     ax2.Visible = 'off';
%     ax2.XTick = [];
%     ax2.YTick = [];
%     
%     % Color 
%     colormap(ax1,'jet')
%     colormap(ax2)
%     
%     % get everthin lined up
%     cb1 = colorbar(ax1,'Position',[0.2 0.1 0.05 0.815]); % four-elements vector to specify Position [left bottom width height]
%     cb2 = colorbar(ax2,'Position',[0.81 0.1 0.05 0.815]);
%     cb1.Label.String = 'Door displacement';
%     cb2.Label.String = 'Cavity acoustic pressure';
%     cb1.Label.FontSize = 14;
%     cb2.Label.FontSize = 14;
%     
% %     view([-70.5 28.505859375]);
% end