close all; clc; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 1 - Calcule as 5 primeiras frequ^encias naturais e formas modais da
% placa e da cavidade quando desacopladas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Part 1 - Plate mode shape
fprintf('Solving plate model\n')
plate_parameters = containers.Map;


%% Part 1 - Cavity mode shape
fprintf('Solving acoustic cavity model\n')
cvt_param = containers.Map;

cvt_param('W') = 240e-3; % volume width [m] (X axis
cvt_param('D') = 240e-3; % volume depth [m] (Y axis)
cvt_param('H') = 540e-3; % volume height [m] (Z axis)
cvt_param('c') = 348; % sound speed [m/s]
cvt_param('rho') = 1.7; % density [kg/m**3]
cvt_param('eSize') = 0.02; % element size [m]

cavity = AcousticModel(cvt_param);
cavity = cavity.get_element_matrices(cavity);
cavity = cavity.generate_mesh(cavity);
cavity = cavity.get_global_matrices(cavity);
spy(cavity.Qg)
