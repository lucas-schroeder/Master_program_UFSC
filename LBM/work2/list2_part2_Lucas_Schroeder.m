%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This is a simple implementation of the LBGK model.
%%% By Andrey R. da Silva, August 2010
%%% Edited by Lucas Schroeder, February 2021
%%% The code does not take into account any specific boundary condition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, clc, close all

fprintf("Muffler\n")
start_time = now;
fprintf("Start Time: "+datestr(start_time, 'yyyy-mm-dd HH:MM:SS')+"\n");

% Block 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lattice size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


r = ceil(30/8);         % resolution [cells/mm] 

Nr = 53*r+2;            % Number of lines   (cells in the y direction)
Mc = 197*r+2;           % Number of columns (cells in the x direction)


% Block 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Physical parameters (macro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_p = 340;                    % Sound velocity on the fluid [m/s]
rho_p = 1.2;                  % physical density [kg/m^3]
nu_STP_p = 1.5e-5;            % Air STP Kinematic Viscosity [m^2/s]
Lx = 0.197;                   % Maximun dimenssion in the x direction [m]
Ly = 0.053;                   % Maximun dimenssion on th y direction  [m]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% definition of the solid static boundaries by vetors using "crossing"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Muffler walls
x_keypoints = r.*[0 0 91 91 128 128 197 197 128 128 91 91 0] + 1;
y_keypoints = r.*[23.35 29.65 29.65 53 53 29.65 29.65 23.35 23.35 0 0 23.35 23.35] + 1;

[vec1,vec2,vec3,vec4,vec5,vec6,vec7,vec8] = crossing3(Nr,Mc,x_keypoints,y_keypoints);


% Block 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lattice parameters (micro - lattice unities)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cs = 1/sqrt(3);                % lattice speed of sound
cs2 = cs^2;                    % Squared speed of sound cl^2
zeta = c_p/cs;                 % Conversion factor between physiscal and lattice units
Dx = Lx/(Mc-2);                % Lattice space (pitch) [m/lattice_length]
Dt = Dx/zeta;                  % time step [s/lattice_time]

omega = 1.9;                   % Relaxation frequency
tau = 1/omega;                 % Relaxation time
nu = cs2*(tau - 0.5);          % Kinematic Viscosity [lattice units]
nu_p = zeta*Dx*nu;             % Physical Kinematic Viscosity [m^2/s]
nu_STP = nu_STP_p/(zeta*Dx);   % STP  Kinematic Viscosity, lattice units
tau_STP = 1/2 + nu_STP/cs2;    % STP Relaxation time
omega_STP = 1/tau_STP;
rho_l = 1;                     % avereged fluid density (latice density)

% Block 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lattice properties for the D2Q9 model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_c = 9 ;                            % number of directions of the D2Q9 model
C_x = [1 0 -1  0 1 -1 -1  1 0];      % velocity vectors in x
C_y = [0 1  0 -1 1  1 -1 -1 0];      % velocity vectors in y
w0 = 4/9. ; w1 = 1/9. ; w2 = 1/36.;  % lattice weights
W  = [w1 w1 w1 w1 w2 w2 w2 w2 w0];   % slide 16
f1 = 3.;
f2 = 4.5;
f3 = 1.5;                            % coef. of the f equil.


% Array of distribution and relaxation functions
f = zeros(Nr,Mc,N_c);                                 
feq = zeros(Nr,Mc,N_c);

% Filling the initial distribution function (at t=0) with initial values
f(:,:,:) = rho_l/9;   
ux = zeros(Nr, Mc);
uy = zeros(Nr, Mc);
% rho_l = 0.01;   % initial disturbance




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ABC parameters (non-reflecting boundary)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = 30;                                     % thickness of the ABC layer
sigma_m = 0.3;                              % constant - empirical
f_t = zeros(Nr,Mc,9);                       % Constructing the target dist. functions
sigr  = zeros(Nr,Mc,9);
sigl  = zeros(Nr,Mc,9);
fs = zeros(Nr,Mc,9);

% targets for the ABC condition
ux_t = 0;
uy_t = 0;
rho_t = 1;

% source distribution functions (uses the Maxweall formula to fill f_t
% based on target values fo rho and u
uxsq = ux_t.^2; 
uysq = uy_t.^2; 
usq = uxsq + uysq;
f_t(:,:,1)= rho_t*w1 .*(1 +f1*ux_t +f2.*uxsq -f3*usq);
f_t(:,:,2)= rho_t*w1 .*(1 +f1*uy_t +f2*uysq -f3*usq);
f_t(:,:,3)= rho_t*w1 .*(1 -f1*ux_t +f2*uxsq -f3*usq);
f_t(:,:,4)= rho_t*w1 .*(1 -f1*uy_t +f2*uysq -f3*usq);
f_t(:,:,5)= rho_t*w2 .*(1 +f1*(+ux_t+uy_t) +f2*(+ux_t+uy_t).^2 -f3.*usq);
f_t(:,:,6)= rho_t*w2 .*(1 +f1*(-ux_t+uy_t) +f2*(-ux_t+uy_t).^2 -f3.*usq);
f_t(:,:,7)= rho_t*w2 .*(1 +f1*(-ux_t-uy_t) +f2*(-ux_t-uy_t).^2 -f3.*usq);
f_t(:,:,8)= rho_t*w2 .*(1 +f1*(+ux_t-uy_t) +f2*(+ux_t-uy_t).^2 -f3.*usq);
f_t(:,:,9)= rho_t*w0 .*(1 - f3*usq);

% OUTLET - right side of the domain
for i=1:9                                   % reconstructing f_t with feq. Notice t    
    for et = 1:Nr    
    sigr(et,Mc-D:Mc,i) = sigma_m*(linspace(0,D,D+1)/(D)).^2;
    end
end

% INLET - left side of the domain
for i=1:9                                   % reconstructing f_t with feq. Notice t    
    for et = 95:119 %1:Nr    
    sigl(et,1:D+1,i) = fliplr(sigma_m*(linspace(0,D,D+1)/(D)).^2);
    end
end

%Block 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin the iteractive process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

total_time = ceil(10*Mc/cs);
time = 0:1:total_time-1;
p_inlet_t = zeros(1,total_time);
ux_inlet_t = zeros(1,total_time);
p_outlet_t = zeros(1,total_time);
ux_outlet_t = zeros(1,total_time);

% cossine sweep
sweep_time = 1:total_time;          % timesteps
f_start = eps*(Dx/zeta);            % frequency at initial time
t_finish = ceil(2*Mc/cs);           % reference time
f_finish = 16e3*(Dx/zeta);          % frequency at reference time
sweep = chirp(sweep_time, f_start, t_finish, f_finish);
sweep(t_finish:total_time)=0;
% plot(sweep_time,sweep);
% title("Chirp excitation");
% xlabel("Time, timesteps");
% ylabel("Density fluctuation, lattice units");

% targets for the inlet boundary condition
ux_inlet = 0;
uy_inlet = 0;
rho_inlet = 1;
feq_inlet = zeros(Nr,Mc,N_c);
uxsq_inlet = ux_inlet.^2; 
uysq_inlet = uy_inlet.^2; 
usq_inlet = uxsq_inlet + uysq_inlet;  

% progress bar
prog = waitbar(0,sprintf('Progress: %0.0f %%',0),...
                      'Name','Running Simulation',...
                      'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(prog,'canceling',0);
for ta = 1 : total_time
    % Check for clicked Cancel button
    if getappdata(prog,'canceling')
        break
    end
%     waitbar(ta/total_time,prog,sprintf('Progress: %0.0f %%',ta/total_time*100));
    waitbar(ta/total_time,prog,sprintf('Timestep %0.0f of %0.0f',ta,total_time));

    % Block 5.1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % propagation (streaming)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    f(:,:,1) = [f(:,1:2,1) f(:,2:Mc-1,1)];
    f(:,:,2) = [f(1:2,:,2);f(2:Nr-1,:,2)];
    f(:,:,3) = [f(:,2:Mc-1,3) f(:,Mc-1:Mc,3)];
    f(:,:,4) = [f(2:Nr-1,:,4);f(Nr-1:Nr,:,4)];
    f(:,:,5) = [f(:,1:2,5) f(:,2:Mc-1,5)];
    f(:,:,5) = [f(1:2,:,5);f(2:Nr-1,:,5)];     % comp ertical
    f(:,:,6) = [f(:,2:Mc-1,6) f(:,Mc-1:Mc,6)];
    f(:,:,6) = [f(1:2,:,6);f(2:Nr-1,:,6)];     % comp fertical
    f(:,:,7) = [f(:,2:Mc-1,7) f(:,Mc-1:Mc,7)];
    f(:,:,7) = [f(2:Nr-1,:,7);f(Nr-1:Nr,:,7)]; % comp ertical
    f(:,:,8) = [f(:,1:2,8) f(:,2:Mc-1,8)];
    f(:,:,8) = [f(2:Nr-1,:,8);f(Nr-1:Nr,:,8)]; % comp ertical
    
    
 
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % restoring the cells that have crossed the static boundary (bounce back)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    G = f; 
    G(vec1) = f(vec3);
    G(vec2) = f(vec4);
    G(vec3) = f(vec1);
    G(vec4) = f(vec2);
    G(vec5) = f(vec7);
    G(vec6) = f(vec8);
    G(vec7) = f(vec5);
    G(vec8) = f(vec6);
    f = G; 
    
    % Block 5.2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % recalculating rho and u
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rho=sum(f,3); 
       
   
    % imposicao do grdiente de pressao para barreira aberta

    %     rho(2:Nr-1,2:3) = 1+.01;%*sin(2*pi*freq*(ta-1));
    %     rho(2:Nr-1,Nr-2:Nr-1) = 1-.01;%*sin(2*pi*freq*(ta-1));


    % Determining the velocities according to Eq.() (see slides)    
    ux = (C_x(1).*f(:,:,1)+C_x(2).*f(:,:,2)+C_x(3).*f(:,:,3)+C_x(4).*f(:,:,4)+C_x(5).*f(:,:,5)+C_x(6).*f(:,:,6)+C_x(7).*f(:,:,7)+C_x(8).*f(:,:,8))./rho ;
    uy = (C_y(1).*f(:,:,1)+C_y(2).*f(:,:,2)+C_y(3).*f(:,:,3)+C_y(4).*f(:,:,4)+C_y(5).*f(:,:,5)+C_y(6).*f(:,:,6)+C_y(7).*f(:,:,7)+C_y(8).*f(:,:,8))./rho ;
    
    % keeping the results
    p_outlet_t(ta) = mean(rho(95:119,Mc-8*r)-rho_l, 'all')*cs2; % isothermal model [lattice uinits]
    ux_outlet_t(ta) = mean(ux(95:119,Mc-8*r), 'all');

    p_inlet_t(ta) = mean(rho(95:119,30*r)-rho_l, 'all')*cs2; % isothermal model [lattice uinits]
    ux_inlet_t(ta) = mean(ux(95:119,30*r), 'all');
    
    % Block 5.3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determining the relaxation functions for each direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    rt0= w0*rho;
    rt1= w1*rho;
    rt2= w2*rho;
    
    uxsq = ux.^2; 
    uysq = uy.^2; 
    usq = uxsq + uysq; 
       
    feq(:,:,1)= rt1 .*(1 +f1*ux +f2.*uxsq -f3*usq);
    feq(:,:,2)= rt1 .*(1 +f1*uy +f2*uysq -f3*usq);
    feq(:,:,3)= rt1 .*(1 -f1*ux +f2*uxsq -f3*usq);
    feq(:,:,4)= rt1 .*(1 -f1*uy +f2*uysq -f3*usq);
    feq(:,:,5)= rt2 .*(1 +f1*(+ux+uy) +f2*(+ux+uy).^2 -f3.*usq);
    feq(:,:,6)= rt2 .*(1 +f1*(-ux+uy) +f2*(-ux+uy).^2 -f3.*usq);
    feq(:,:,7)= rt2 .*(1 +f1*(-ux-uy) +f2*(-ux-uy).^2 -f3.*usq);
    feq(:,:,8)= rt2 .*(1 +f1*(+ux-uy) +f2*(+ux-uy).^2 -f3.*usq);
    feq(:,:,9)= rt0 .*(1 - f3*usq);
    
    % Implement Bondary Condition - Cossine sweep at entrance, ABC
    rho_inlet = 1 + 0.1*sweep(ta);
    
    rt0_inlet= w0*rho_inlet;
    rt1_inlet= w1*rho_inlet;
    rt2_inlet= w2*rho_inlet;
    
    feq_inlet(:,:,1)= rt1_inlet.*(1 +f1*ux_inlet +f2.*uxsq_inlet -f3*usq_inlet);
    feq_inlet(:,:,2)= rt1_inlet.*(1 +f1*uy_inlet +f2*uysq_inlet -f3*usq_inlet);
    feq_inlet(:,:,3)= rt1_inlet.*(1 -f1*ux_inlet +f2*uxsq_inlet -f3*usq_inlet);
    feq_inlet(:,:,4)= rt1_inlet.*(1 -f1*uy_inlet +f2*uysq_inlet -f3*usq_inlet);
    feq_inlet(:,:,5)= rt2_inlet.*(1 +f1*(+ux_inlet+uy_inlet) +f2*(+ux_inlet+uy_inlet).^2 -f3.*usq_inlet);
    feq_inlet(:,:,6)= rt2_inlet.*(1 +f1*(-ux_inlet+uy_inlet) +f2*(-ux_inlet+uy_inlet).^2 -f3.*usq_inlet);
    feq_inlet(:,:,7)= rt2_inlet.*(1 +f1*(-ux_inlet-uy_inlet) +f2*(-ux_inlet-uy_inlet).^2 -f3.*usq_inlet);
    feq_inlet(:,:,8)= rt2_inlet.*(1 +f1*(+ux_inlet-uy_inlet) +f2*(+ux_inlet-uy_inlet).^2 -f3.*usq_inlet);
    feq_inlet(:,:,9)= rt0_inlet.*(1 - f3*usq_inlet);
    
    
   
    % Block 5.4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Collision (relaxation) step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    f = (1-omega)*f + omega*feq - sigr.*(feq - f_t) - sigl.*(feq-feq_inlet);

    
 
%     % Ploting the results in real time  \
%     figure(1)
%     surf(rho-1), view(2), shading flat, axis equal, caxis([-.001 .001])
%     % surf(rho-1), view(2), shading flat, axis equal, caxis([-.001 .001])
%     % colormap(gray)
%     grid on
%     % figure(2)
%     % plot(ux(150,:))
%     % figure(3)
%     % surf(sqrt(ux.^2+uy.^2)), view(2), shading flat, axis equal, caxis([-.01 .01])
% %     pause(.0001)    

end %  End main time Evolution Loop
delete(prog); % close waitbar
end_time = now;
fprintf("End time: "+datestr(end_time, 'yyyy-mm-dd HH:MM:SS')+"\n");
fprintf("Simulation Elapsed Time: "+datestr(end_time-start_time, 'HH:MM:SS')+"\n");
writematrix([time' p_inlet_t' p_outlet_t'],"pressure_history"+datestr(end_time, 'yyyy-mm-dd-HHMMSS')+".csv");
%% Plot results
figure(1)
plot(time,p_inlet_t);hold on;
plot(time,p_outlet_t);hold off;
title("Inlet and Outlet pressure");
xlabel("Time, timesteps");
ylabel("Pressure, lattice units");
legend("p_{inlet}","p_{outlet}");
% xlim([0 1350]);
writematrix([time' p_inlet_t' p_outlet_t'],"pressure_history"+datestr(now, 'yyyy-mm-dd-HHMMSS')+".csv");

figure(2)
% surf(rho-1), view(2), shading flat, axis equal, caxis([-0.2 0.2]), colorbar;

surf(ux), view(2), shading flat, axis equal, caxis([-0.001 0.001]), colorbar;
hold on;
plot(Mc-8*r,95:119, "xk");
plot(30*r,95:119, "xk");
plot(x_keypoints,y_keypoints,"k");
hold off;
colormap jet;
title("Last timestep u_x countour plot");
xlim([0 790]);
ylim([0 214]);

%% Calculate Transmission Loss
Pin_f = fft(p_inlet_t);
Pout_f = fft(p_outlet_t);
L = length(p_outlet_t);
H = 20*log10(abs(Pout_f./Pin_f));

Df = 1/Dt;                 % sampling frequency, Hz
time = Dt*(0:(L - 1));     % time vector, s
frq = Df*(0:(L - 1)/2)/L;  % frequency vector, Hz

results_FEM = openfig('CurvaMufflerMareze.fig'); % load and plot FEM results
hold on;
plot(frq, H(1:length(frq)),"Xr");
title("Transfer function");
ylabel("H(f), dB");
xlabel("Frequency, Hz");
legend("FEM","LBGK");
grid on;
xlim([0 10e3]);
ylim([-50 20]);



