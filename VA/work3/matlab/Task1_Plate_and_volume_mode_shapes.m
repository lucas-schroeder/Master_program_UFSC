close all; clc; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 1 - Calcule as 5 primeiras frequ^encias naturais e formas modais da
% placa e da cavidade quando desacopladas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Part 1 - Plate mode shape
fprintf('Solving plate model\n')
plate_param = containers.Map;
plate_param('h') = 1.6e-3;  % thickness [m]
plate_param('H') = 500e-3;  % door height [m]
plate_param('L') = 200e-3;  % door length [m]
plate_param('E') = 186e9;  % modulus of elasticity [Pa]
plate_param('nu') = 0.3;  % poison coef
plate_param('rho') = 7870; % steel density [kg/m**3]
plate_param('eSize') = 0.02; % element size [m]

plate = PlateModel(plate_param);
plate = plate.get_element_matrices(plate);
plate = plate.generate_mesh(plate);
plate = plate.get_global_matrices(plate);
spy(plate.M)

% Two-points recatangle
regions=containers.Map;
regions('1')=[0.000 0.100, 0.200 0.120];
regions('2')=[0.000 0.380, 0.200 0.400];
regions('3')=[0.180 0.260, 0.220 0.280];

plate = plate.solve_eigenvalue_problem(plate, 50);

%% Part 2 - Cavity mode shape
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
Vc = cvt.results.Phi;

%% 2.1 - Plot decoupled acoustic cavity mode shapes
xp1 = 0:cvt_param('eSize'):cvt_param('W');
yp1 = 0:cvt_param('eSize'):cvt_param('D');
zp1 = 0:cvt_param('eSize'):cvt_param('H');

for mode = 1:5
    [X1,Y1,Z1] = meshgrid(xp1,yp1,zp1);
    figure()
    set(gcf, 'Position', get(0, 'Screensize'));
    mode_shape = reshape(Vc(:,mode), [cvt.nNodesX,cvt.nNodesY,cvt.nNodesZ]);
    mode_shape = permute(mode_shape,[2 1 3]);
    slice(X1,Y1,Z1, mode_shape, 0, 0, 0.540);
    box on
    ax = gca;
    ax.BoxStyle = 'full';
    shading interp
    set(gca, 'FontSize', 20)
    xlabel('Comprimento [m]')
    ylabel('Largura [m]')
    zlabel('Altura [m]')
    title(sprintf('Decoupled acoustic cavity - Mode %.f (%0.1f [Hz])',mode,cvt.results.fn(mode)));
    colormap jet;colorbar
    set(gcf,'color','w');
    axis equal
end
