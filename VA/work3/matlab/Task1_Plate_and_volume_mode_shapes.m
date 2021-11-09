close all; clc; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Work 3 - Computational Vibroacoustics
% Lucas Schroeder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
%% Part 1 - Plate mode shape %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('1 - Solving door model\n')
door_param = containers.Map;
door_param('h') = 2e-3;  % Thickness [m]
door_param('H') = 500e-3;  % Door height [m]
door_param('L') = 200e-3;  % Door length [m]
door_param('E') = 200e9;   % Modulus of elasticity [Pa]
door_param('nu') = 0.3;    % Poison coef
door_param('rho') = 7850;  % Steel density [kg/m**3]
door_param('eSize') = 0.02; % Element size [m]

door = PlateModel(door_param);
door = door.get_element_matrices(door);
door = door.generate_mesh(door);
door = door.get_global_matrices(door);

% Fixed regions (Two-points recatangles [x1 z1, x2 z2])
regions = containers.Map;
regions('1') = [0.000 0.100, 0.020 0.120]; % Lower hinge
regions('2') = [0.000 0.380, 0.020 0.400]; % Upper hinge
regions('3') = [0.180 0.240, 0.200 0.280]; % Knob

door = door.apply_bc(door, regions);
door = door.solve_eigenvalue_problem(door, 50);
fprintf('Done!\n')

%% 1.1 - Plot door mode shapes (decoupled)
x = 0 : door.dx : door.L;
z = 0 : door.dz : door.H;
[X, Z] = meshgrid(x, z);
figure()
for mode=1:6
    subplot(2, 3, mode)
    mode_shape = door.results.Vc(:, mode);    
    mode_shape = 0.02.*mode_shape./max(mode_shape); %scale
    mode_shape = reshape(mode_shape, door.nNodesX, door.nNodesZ);
    surf(X, mode_shape.', Z, mode_shape.');
    title(sprintf('Mode %.f (%0.1f [Hz])',mode,door.results.fn(mode)))
    xlabel('X [m]');
    zlabel('Z [m]');
    set(gca,'ytick',[])
    axis equal;
    grid on;
    box on;
end


%% Part 2 - Solve acoutic cavity model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('2 - Solving acoustic cavity model\n')

cvt_param = containers.Map;
cvt_param('W') = 240e-3; % volume width [m] (X axis
cvt_param('D') = 240e-3; % volume depth [m] (Y axis)
cvt_param('H') = 540e-3; % volume height [m] (Z axis)
cvt_param('c') = 346.25; % sound speed [m/s]
cvt_param('rho') = 1.225; % density [kg/m**3]
cvt_param('eSize') = 0.02; % element size [m]

cvt = AcousticModel(cvt_param);
cvt = cvt.get_element_matrices(cvt);
cvt = cvt.generate_mesh(cvt);
cvt = cvt.get_global_matrices(cvt);
cvt = cvt.solve_eigenvalue_problem(cvt, 25);
fprintf('Done!\n')


%% 2.1 - Plot acoustic cavity mode shapes (decoupled)
Vc_cpl = cvt.results.Phi;
x = 0 : cvt.dx : cvt.W;
y = 0 : cvt.dy : cvt.D;
z = 0 : cvt.dz : cvt.H;
[X, Y, Z] = meshgrid(x, y, z);
xslice = [0.000, cvt.W];
yslice = [0.000, cvt.D];
zslice = [0.000, cvt.H];

figure()
for mode = 1:6
    subplot(2, 3, mode)
    mode_shape = reshape(Vc_cpl(:,mode), [cvt.nNodesX, cvt.nNodesY, cvt.nNodesZ]);
    mode_shape = permute(mode_shape,[2 1 3]);
    slice(X,Y,Z, mode_shape, xslice, yslice, zslice);
    box on
    ax = gca;
    ax.BoxStyle = 'full';
    shading interp
    xlabel('X [m]')
    ylabel('Y [m]')
    zlabel('Z [m]')
    title(sprintf('Mode %.f (%0.1f [Hz])',mode,cvt.results.fn(mode)));
    colormap jet;
    set(gcf,'color','w');
    axis equal
end
c = colorbar;
c.Label.String = 'Acoustic Pressure';
c.Position = [0.9 0.1 0.03 0.8];

%% Part 3 - Solve acoutic cavity with back opening %%%%%%%%%%%%%%%%%%%%%%%%
fprintf('3 - Solve acoutic cavity with back opening\n')

cvt2 = AcousticModel(cvt_param);
cvt2 = cvt2.get_element_matrices(cvt2);
cvt2 = cvt2.generate_mesh(cvt2);
cvt2 = cvt2.get_global_matrices(cvt2);

% Constant pressure regions (Two-points box [x1 y1 z1, x2 y2 z2])
regions = containers.Map;
regions('1')=[0.180 0.240 0.040; 0.220 0.240 0.100]; % back opening

cvt2 = cvt2.apply_bc(cvt2, regions);
cvt2 = cvt2.solve_eigenvalue_problem(cvt2, 25);
fprintf('Done!\n')

%% 3.1 - Plot acoutic cavity with opening
Vc_open = cvt2.results.Phi;
x = 0 : cvt.dx : cvt.W;
y = 0 : cvt.dy : cvt.D;
z = 0 : cvt.dz : cvt.H;
[X, Y, Z] = meshgrid(x, y, z);
xslice = [0.000, cvt.W];
yslice = [0.000, cvt.D];
zslice = [0.000, cvt.H];

figure()
for mode = 1:6    
    subplot(2,3,mode)
    mode_shape = reshape(Vc_open(:,mode), [cvt.nNodesX, cvt.nNodesY, cvt.nNodesZ]);
    mode_shape = permute(mode_shape,[2 1 3]);
    slice(X,Y,Z, mode_shape, xslice, yslice, zslice);
    box on;
    grid on;
    ax = gca;
    ax.BoxStyle = 'full';
    shading interp
    xlabel('X [m]')
    ylabel('Y [m]')
    zlabel('Z [m]')
    title(sprintf('Mode %.f (%0.1f [Hz])',mode,cvt2.results.fn(mode)));
    colormap jet;
    set(gcf,'color','w');
    axis equal
end
c = colorbar;
c.Label.String = 'Acoustic Pressure';
c.Position = [0.9 0.1 0.03 0.8];

%% Part 4 - Coupling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('4 - Calculating coupling matrix\n')
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
% The door and the cavity meshes are not coincident, so we must determine 
% their relation. Using cavity mesh coordinate system (csys) as the default
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
        fprintf('Meshes are not concident!!\n')
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
        fprintf('Meshes are not concident!!\n')
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
            % door_dof(z) is the i-th door DoF in the door csys
            door_dof(z) = door.node_dof(J, i);
            z = z + 1;
        end          
    end
    
    r = 1; 
    e2 = door_cvt_elements(e); % Convert the element from door csys to cavity csys
    for j = 1:8 % For each node of the Hex element
        J = cvt.element_con(e2, j);  % Global node (in cavity csys) 
        for i = 1
            % cvt_dof(z) is the i-th cavity DoF in the cavity csys
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
mask = true(door.ndof, 1); % Boolean mask for the rows
mask(door.fixed_dof) = false;
S = S(mask, :);
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
fprintf('...Solving the coupled system eigenvalue problem\n');
opts.v0 = (1:length(K_coupled)).'*200; % Initial state
[Vc_cpl, Wn2] = eigs(K_coupled, M_coupled, 50, 'sm', opts);
fn_cpl = diag(Wn2).^(0.5)./(2*pi); % in Hertz

% Normalizing eigenvectors matrix by the mass matrix,
% such that Vc' * M * Vc = I
% m_r = abs(diag(Vc_cpl'*M_coupled*Vc_cpl));
% m_r = 1./(m_r.^0.5);
% for i = 1:size(Vc_cpl,2)
%     Vc_cpl(:,i) = Vc_cpl(:,i).*m_r(i);
% end

% Ignore complex residues
fn_cpl = real(fn_cpl);
Vc_cpl = abs(Vc_cpl);

% Ordering eigenvalues and eigenvectors
[fn_cpl, idx] = sort(fn_cpl);
Vc_cpl = Vc_cpl(:, idx);
fprintf('Done!\n');


%% 4.2 - Splitting the DoF vector
Vc_door = Vc_cpl(1:length(door.M), :);
% Normalizing eigenvectors matrix by the mass matrix,
% such that Vc' * M * Vc = I
m_r = abs(diag(Vc_door'*door.M*Vc_door));
m_r = 1./(m_r.^0.5);
for i = 1:size(Vc_door,2)
    Vc_door(:,i) = Vc_door(:,i).*m_r(i);
end

% Adding fixed Dof back into the displacement vector
nModos = size(Vc_door, 2);
fixed = zeros(1, nModos);
for row = door.fixed_dof
    Vc_door = [Vc_door(1:row-1,:) ; fixed ; Vc_door(row:end,:)];
end

% Taking only the vertical displacements DoF from door
Vc_door = real(Vc_door(1:3:end,:)); % Displacement shape (dof 1,4,7,10,13,...)

% Taking the cavity results
Vc_cvt = Vc_cpl(length(door.M)+1:end, :);
% Normalizing eigenvectors matrix by the mass matrix,
% such that Vc' * M * Vc = I
m_r = abs(diag(Vc_cvt'*cvt.Qg*Vc_cvt));
m_r = 1./(m_r.^0.5);
for i = 1:size(Vc_cvt,2)
    Vc_cvt(:,i) = Vc_cvt(:,i).*m_r(i);
end

%% 4.3 - Plot door mode shapes (coupled)
x = 0 : door.dx : door.L;
z = 0 : door.dz : door.H;
[X, Z] = meshgrid(x, z);
figure()
for mode = 1:12
    subplot(3,4,mode)
    mode_shape = Vc_door(:,mode);
    mode_shape = 0.03.*mode_shape./max(mode_shape); %scale
    mode_shape = reshape(mode_shape, door.nNodesX, door.nNodesZ);
    surf(X, mode_shape.', Z, mode_shape.');
    title(sprintf('Mode %.f (%0.1f [Hz])', mode, fn_cpl(mode+1)))
    xlabel('X [m]');
    zlabel('Z [m]');
    set(gca,'ytick',[])
    axis equal;
    grid on;
    box on;
end

%% 4.4 - Plot acoustic cavity mode shapes (coupled)
x = 0 : cvt.dx : cvt.W;
y = 0 : cvt.dy : cvt.D;
z = 0 : cvt.dz : cvt.H;
[X, Y, Z] = meshgrid(x, y, z);
xslice = [0.000, cvt.W];
yslice = [0.000, cvt.D];
zslice = [0.000, cvt.H];

figure()
for mode = 1:12    
    subplot(3,4,mode)
    mode_shape = reshape(Vc_cvt(:,mode), [cvt.nNodesX, cvt.nNodesY, cvt.nNodesZ]);
    mode_shape = permute(mode_shape,[2 1 3]);
    slice(X,Y,Z, mode_shape, xslice, yslice, zslice);
    box on;
    grid on;
    ax = gca;
    ax.BoxStyle = 'full';
    shading interp
    xlabel('X [m]')
    ylabel('Y [m]')
    zlabel('Z [m]')
    title(sprintf('Mode %.f (%0.1f [Hz])', mode, fn_cpl(mode)));
    colormap jet;
    set(gcf,'color','w');
    axis equal
end
c = colorbar;
c.Label.String = 'Acoustic Pressure';
c.Position = [0.9 0.1 0.03 0.8];

%% 4.4b - Plot MAC
% Door MAC
Phi1 = door.results.Vc(:,1:10); % column vector
Phi2 = Vc_door(:,1:10);
MAC_door = zeros(size(Phi1,2),size(Phi2,2));
for j = 1:size(Phi1,2)
    for k = 1:size(Phi2,2)
        MAC_door(j,k) = abs(Phi1(:,j).'*Phi2(:,k)).^2/...
                       ((Phi1(:,j).'*Phi1(:,j))*(Phi2(:,k).'*Phi2(:,k)));
    end
end

figure
imagesc(MAC_door)
title('Modal Assurance Criteria')
xlabel('Decoupled Plate Mode')
ylabel('Coupled Plate Mode')
colorbar

%% 4.5 - Find nodes and DoF of interest for FRF
A_coor = [0.160 0.420]; % Force at the door. Door csys [x z]
B_coor = [0.100 0.160 0.340]; % Response inside cavity. Cavity csys [x y z]

%%%% 
x = door.node_coor(:,1) == A_coor(1);
z = door.node_coor(:,2) == A_coor(2);
idx = x & z;
if sum(idx)==1
    A_node = find(idx==1);
else
    fprintf('There is no node on the selected location\n')
end
A_dof = door.node_dof(A_node,2); % In door csys (also in coupled-system csys)

%%%%
x = cvt.node_coor(:,1) == B_coor(1);
y = cvt.node_coor(:,2) == B_coor(2);
z = cvt.node_coor(:,3) == B_coor(3);
idx = x & y & z;
if sum(idx)==1
    B_node = find(idx==1);
else
    fprintf('There is no node on the selected location\n')
end
B_dof = cvt.node_dof(B_node,2); % In cavity csys
B_dof_coupled_csys = B_dof + length(M); % In coupled-system csys

%% 4.6 - FRF - Weak coupling
fprintf('4.6 - Computing weak coupling FRF\n')
% FRF between a force at poit A a displacements at every door node
freq = 10:2:700; % [Hz]
H_jk = zeros(door.ndof/3, length(freq)); % [m/N]
V_jk = zeros(door.ndof/3, length(freq)); % [(m/s)/N] 
eta = 0.03;
k = find(door.node_dof(:,2) == A_dof); % DoF index without rotations
k = int16(k);
for j = 1:size(door.results.Vc, 1) % For each response dof
    for p = 1:length(freq) % At a given frequency
        for n = 1:length(door.results.fn) % Sum the contribution of every mode
            H_jk(j, p) = H_jk(j, p) + door.results.Vc(j, n)*door.results.Vc(k, n)/...
                                    (door.results.fn(n)^2 - freq(p)^2 + 1j*eta*door.results.fn(n)^2);
            V_jk(j, p) = 1j*freq(p)*H_jk(j, p);
        end
    end
end


% FRF (Pa/(m^3/s/m^3)) between a volumetric acoustic sources at every cavity
% node in contact with theplate, and point B
G_jk = zeros(cvt.ndof, length(freq));
k = find(cvt.node_dof(:,2) == B_dof); % DoF index without rotations
k = int16(k);
for j = 1:size(cvt.results.Phi, 1) % For each response dof
    for p = 1:length(freq) % At a given frequency
        for n = 1:length(cvt.results.fn) % Sum the contribution of every mode
            G_jk(j, p) = G_jk(j, p) + cvt.results.Phi(j, n)*cvt.results.Phi(k, n)/...
                                          (cvt.results.fn(n)^2 - freq(p)^2);
        end
    end
end

G_jk_door = G_jk(door_cvt_nodes, :); % FRF between point B and the DoFs at the door (Pa/(m^3/m^3))
vol_element = cvt.dx*cvt.dy*cvt.dz;
srf_element = cvt.dx*cvt.dz;
P_jk_door = G_jk_door*vol_element/srf_element; % Pa/m/s (Pa per unit displacement rate of element face.

I = zeros(size(H_jk));
for j = 1:size(P_jk_door,1)
    I(j,:) = V_jk(j,:) .* P_jk_door(j,:); % [(m/s)/N]*[Pa/(m/s)] = [Pa/N]
end

total_FRF = sum(real(I), 1);
p_ref = 20*e-6; % 20 uPa
total_FRF_db = 20.*log(abs(total_FRF/p_ref));
fprintf('Done!\n')

%% 4.7 - Plot ansys and code
fprintf('4.7 - Importing Ansys results for comparison\n')
Ansys_pointA = xlsread('Ansys_pointA.xlsx'); % Meter/Newton
Ansys_pointA = [Ansys_pointA(:,2) Ansys_pointA(:,3).*exp(1j*deg2rad(Ansys_pointA(:,4)))];
Ansys_pointB = xlsread('Ansys_pointB.xlsx'); % Pascal(B)/Newton(A)
Ansys_pointB = [Ansys_pointB(:,2) Ansys_pointB(:,3).*exp(1j*deg2rad(Ansys_pointB(:,4)))];

figure(); hold on;
f = real(Ansys_pointB(:,1));
plot(f, 20.*log(abs(Ansys_pointB(:,2))/p_ref))
plot(freq, total_FRF_db)
legend({'Ansys', 'Author'})
xlabel('Frequency [Hz]') 
ylabel('SPL [dB] (ref. 20 \muPa)')
fprintf('Done!\n')

%% 4.8 - FRF - Strong coupling
fprintf('4.8 - Computing strong coupling FRF\n')
% FRF between a force at poit A a displacements at every door node
freq = 10:2:700; % [Hz]
H_jk_str = zeros(door.ndof/3, length(freq)); % [m/N]
V_jk_str = zeros(door.ndof/3, length(freq)); % [(m/s)/N] 
eta = 0.03;
k = find(door.node_dof(:,2) == A_dof); % DoF index without rotations
k = int16(k);
for j = 1:size(Vc_door, 1) % For each response dof
    for p = 1:length(freq) % At a given frequency
        for n = 1:length(fn_cpl) % Sum the contribution of every mode
            H_jk_str(j, p) = H_jk_str(j, p) + Vc_door(j, n)*Vc_door(k, n)/...
                                    (fn_cpl(n)^2 - freq(p)^2 + 1j*eta*fn_cpl(n)^2);
            V_jk_str(j, p) = 1j*freq(p)*H_jk_str(j, p);
        end
    end
end


% FRF (Pa/(m^3/s/m^3)) between a volumetric acoustic sources at every cavity
% node in contact with theplate, and point B
G_jk_str = zeros(cvt.ndof, length(freq));
k = find(cvt.node_dof(:,2) == B_dof); % DoF index without rotations
k = int16(k);
for j = 1:size(Vc_cvt, 1) % For each response dof
    for p = 1:length(freq) % At a given frequency
        for n = 1:length(fn_cpl) % Sum the contribution of every mode
            G_jk_str(j, p) = G_jk_str(j, p) + Vc_cvt(j, n)*Vc_cvt(k, n)/...
                                    (fn_cpl(n)^2 - freq(p)^2);
        end
    end
end

G_jk_door_str = G_jk_str(door_cvt_nodes, :); % FRF between point B and the DoFs at the door (Pa/(m^3/m^3))
vol_element = cvt.dx*cvt.dy*cvt.dz;
srf_element = cvt.dx*cvt.dz;
P_jk_door_str = G_jk_door_str*vol_element/srf_element; % Pa/m/s (Pa per unit displacement rate of element face.

I_str = zeros(size(H_jk));
for j = 1:size(P_jk_door,1)
    I_str(j,:) = V_jk_str(j,:) .* P_jk_door_str(j,:); % [(m/s)/N]*[Pa/(m/s)] = [Pa/N]
end

total_FRF_strong = sum(real(I), 1);
p_ref = 20*e-6; % 20 uPa
total_FRF_db_str = 20.*log(abs(total_FRF_strong/p_ref));
fprintf('Done!\n')

figure(); hold on;
plot(freq, total_FRF_db_str)
title('Strong coupling')
xlabel('Frequency [Hz]') 
ylabel('SPL [dB] (ref. 20 \muPa)')

%% 4.9 - FRF - Weak coupling - with opening
fprintf('4.9 - Computing weak coupling FRF - with opening\n')
% FRF between a force at poit A a displacements at every door node
freq = 10:2:700; % [Hz]
H_jk_open = zeros(door.ndof/3, length(freq)); % [m/N]
V_jk_open = zeros(door.ndof/3, length(freq)); % [(m/s)/N] 
eta = 0.03;
k = find(door.node_dof(:,2) == A_dof); % DoF index without rotations
k = int16(k);
for j = 1:size(door.results.Vc, 1) % For each response dof
    for p = 1:length(freq) % At a given frequency
        for n = 1:length(door.results.fn) % Sum the contribution of every mode
            H_jk_open(j, p) = H_jk_open(j, p) + door.results.Vc(j, n)*door.results.Vc(k, n)/...
                                    (door.results.fn(n)^2 - freq(p)^2 + 1j*eta*door.results.fn(n)^2);
            V_jk_open(j, p) = 1j*freq(p)*H_jk_open(j, p);
        end
    end
end


% FRF (Pa/(m^3/s/m^3)) between a volumetric acoustic sources at every cavity
% node in contact with theplate, and point B
G_jk_open = zeros(cvt2.ndof, length(freq));
k = find(cvt2.node_dof(:,2) == B_dof); % DoF index without rotations
k = int16(k);
for j = 1:size(cvt2.results.Phi, 1) % For each response dof
    for p = 1:length(freq) % At a given frequency
        for n = 1:length(cvt2.results.fn) % Sum the contribution of every mode
            G_jk_open(j, p) = G_jk_open(j, p) + cvt2.results.Phi(j, n)*cvt2.results.Phi(k, n)/...
                                          (cvt2.results.fn(n)^2 - freq(p)^2);
        end
    end
end

G_jk_door_open = G_jk_open(door_cvt_nodes, :); % FRF between point B and the DoFs at the door (Pa/(m^3/m^3))
vol_element = cvt2.dx*cvt2.dy*cvt2.dz;
srf_element = cvt2.dx*cvt2.dz;
P_jk_door_open = G_jk_door_open*vol_element/srf_element; % Pa/m/s (Pa per unit displacement rate of element face.

I_open = zeros(size(H_jk));
for j = 1:size(P_jk_door_open,1)
    I_open(j,:) = V_jk(j,:) .* P_jk_door_open(j,:); % [(m/s)/N]*[Pa/(m/s)] = [Pa/N]
end

total_FRF_open = sum(real(I_open), 1);
p_ref = 20*e-6; % 20 uPa
total_FRF_db_open = 20.*log(abs(total_FRF_open/p_ref));
fprintf('Done!\n')

figure(); hold on;
plot(freq, total_FRF_db_open)
plot(freq, total_FRF_db)
legend('Open cavity','Closed cavity')
title('Weak coupling - open cavity')
xlabel('Frequency [Hz]') 
ylabel('SPL [dB] (ref. 20 \muPa)')

%% 4.10 - Phase speed of bending waves in thin plates

h = door_param('h'); % Thickness [m]
E = door_param('E'); % Modulus of elasticity [Pa]
nu = door_param('nu'); % Poison coef
rho = door_param('rho'); % Steel density [kg/m**3]
f_min = 10; w_min = 2*pi*f_min;
f_max = 700; w_max = 2*pi*f_max;

c_min = sqrt(w_min*sqrt(E*h^3/(12*rho*(1-nu^2))));
c_max = sqrt(w_max*sqrt(E*h^3/(12*rho*(1-nu^2))));

lambda_min = c_min / f_min;
lambda_max = c_max / f_max;

fprintf('Smallest expected bending wavelength: %0.1f mm \n', 1000*lambda_max)
fprintf('Smallest expected air wavelength: %0.1f mm \n', 1000*cvt_param('c')/f_max)

total_time = toc;
fprintf('Total computation time: %0.1f seconds\n', total_time)