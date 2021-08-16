%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This is a simple implementation of the LBGK model.
%%% By Andrey R. da Silva, August 2010
%%% Edited by Lucas Schroeder, February 2021
%%% The code does not take into account any specific boundary condition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, clc, close all

fprintf("Axisymmetric Zhou et al. Model - Without Flow\n")
start_time = now;
fprintf("Start Time: "+datestr(start_time, 'yyyy-mm-dd HH:MM:SS')+"\n");

% Block 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lattice size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nr = 252;            % Number of lines   (cells in the y direction)
Mc = 502;            % Number of columns (cells in the x direction)
r_a = 20;            % Duct radius


% Block 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Physical parameters (macro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_p = 340;                    % Sound velocity on the fluid [m/s]
rho_p = 1.2;                  % physical density [kg/m^3]
nu_STP_p = 1.5e-5;            % Air STP Kinematic Viscosity [m^2/s]
Lx = 0.500;                   % Maximun dimenssion in the x direction [m]
Ly = 0.250;                   % Maximun dimenssion on th y direction  [m]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% definition of the solid static boundaries by vetors using "crossing"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cylinder Walls
x_cylynder = [31 31 252] + 0.5;
y_cylinder = [0 21 21] + 0.5;
[vec1,vec2,vec3,vec4,vec5,vec6,vec7,vec8] = crossing3(Nr,Mc,x_cylynder,y_cylinder);

% Axis
x_axis = [1 502];
y_axis = [1 1];
[vec12,vec22,vec32,vec42,vec52,vec62,vec72,vec82] = crossing3(Nr,Mc,x_axis,y_axis);

% buffer zones
x_buffer = [31 31 470 470] + 0.5;
y_buffer = [0 220 220 0] + 0.5;

% inlet zones
x_inlet = [61 61] + 0.5;
y_inlet = [0 21] + 0.5;

% Block 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lattice parameters (micro - lattice unities)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cs = 1/sqrt(3);                % lattice speed of sound
cs2 = cs^2;                    % Squared speed of sound cl^2
zeta = c_p/cs;                 % Conversion factor between physiscal and lattice units
Dx = Lx/(Mc-2);                % Lattice space (pitch) [m/lattice_length]
Dt = Dx/zeta;                  % time step [s/lattice_time]

omega = 1.3;                   % Relaxation frequency
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
h_1 = zeros(Nr,Mc,N_c); % Reis et al. (2009) firt order source term
h_2 = zeros(Nr,Mc,N_c); % Reis et al. (2009) second order source term

% Filling the initial distribution function (at t=0) with initial values
f(:,:,:) = rho_l/9;   
ux = zeros(Nr, Mc);
uy = zeros(Nr, Mc);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ABC parameters (non-reflecting boundary)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = 30;                                     % thickness of the ABC layer
sigma_m = 0.3;                              % constant - empirical
f_t = zeros(Nr,Mc,9);                       % Constructing the target dist. functions
sigRight  = zeros(Nr,Mc,9);
sigLeft   = zeros(Nr,Mc,9);
sigTop    = zeros(Nr,Mc,9);
sigInlet  = zeros(Nr,Mc,9);
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

% RIGHT side of the domain
for i=1:9                                   % reconstructing f_t with feq. Notice t    
    for et = 1:Nr    
    sigRight(et,Mc-D:Mc,i) = sigma_m*(linspace(0,D,D+1)/(D)).^2;
    end
end


% LEFT side of the domain
for i=1:9                                       
    for et = 1:Nr    
    sigLeft(et,1:D+1,i) = fliplr(sigma_m*(linspace(0,D,D+1)/(D)).^2);
    end
end


% TOP side of the domain
for i=1:9                                       
    for et = 1:Mc
        sigTop(Nr-D:Nr,et,i) = transpose(sigma_m*(linspace(0,D,D+1)/(D)).^2);
    end
end


% INLET - left side of the duct
for i=1:9                                       
    for et = 2:21
    sigInlet(et,(D+2):(2*D+2),i) = fliplr(sigma_m*(linspace(0,D,D+1)/(D)).^2);
    end
end

% Plot ABS Layers
figure(1)
layer=3;
surf(sigRight(:,:,layer)), view(2), shading flat, colorbar; hold on;
surf(sigLeft(:,:,layer)), view(2), shading flat, colorbar;
surf(sigTop(:,:,layer)), view(2), shading flat, colorbar;
surf(sigInlet(:,:,layer)), view(2), shading flat, colorbar;
plot(x_cylynder,y_cylinder,"k");
plot(x_buffer,y_buffer,"--k");
plot(x_inlet,y_inlet,"--k");
hold off;
title("ABS regions - sigma colormap");
xlim([1 Mc]);
ylim([1 Nr]);
set(gca,"fontsize", 22);
colormap cool;


%Block 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin the iteractive process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_time = ceil(10*250/cs); % 10 wavepass over the duct
transient = 10;
time = 0:1:(total_time-1);


% cossine sweep
sweep_duration = ceil(3/4*(total_time-transient));
sweep_time = 1:sweep_duration;           % timesteps vector
f_start = eps;                           % frequency at initial time
t_finish = sweep_duration;               % reference time
ka_max = 1.8;                            % maximum Helmholtz number
f_finish = cs*ka_max/(2*pi*r_a);         % frequency at the end of sweep
sweep_aux = chirp(sweep_time, f_start, t_finish, f_finish, "linear", -90);
exp_win = [ones(1,length(sweep_aux)-200) exp(-0.02.*(1:1:200))];
sweep = zeros(1,total_time);
sweep(transient:transient+sweep_duration-1) = sweep_aux.*exp_win;

figure(2)
plot(time,0.1*sweep);
    title("Chirp excitation");
    xlabel("Time, timesteps");
    ylabel("\rho ', lattice units");
    ylim([-0.2 0.2]);
    xlim([0 total_time]);
    set(gca,"fontsize", 22);
%%
% targets for the inlet boundary condition
ux_inlet = 0; % FIX ME to include flow
uy_inlet = 0;
rho_inlet = 1;
feq_inlet = zeros(Nr,Mc,N_c);
uxsq_inlet = ux_inlet.^2; 
uysq_inlet = uy_inlet.^2; 
usq_inlet = uxsq_inlet + uysq_inlet;

% Generate measurements vector
rho_outlet_t = zeros(1,total_time);
rho_inlet_t = zeros(1,total_time);
ux_outlet_t= zeros(1,total_time);

% Distance from axis matrix
R = ones(Nr,Mc);
R(1:Nr,:) = R.*[0:251]';

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
    f(:,:,5) = [f(1:2,:,5);f(2:Nr-1,:,5)];     % comp vertical
    f(:,:,6) = [f(:,2:Mc-1,6) f(:,Mc-1:Mc,6)];
    f(:,:,6) = [f(1:2,:,6);f(2:Nr-1,:,6)];     % comp vertical
    f(:,:,7) = [f(:,2:Mc-1,7) f(:,Mc-1:Mc,7)];
    f(:,:,7) = [f(2:Nr-1,:,7);f(Nr-1:Nr,:,7)]; % comp vertical
    f(:,:,8) = [f(:,1:2,8) f(:,2:Mc-1,8)];
    f(:,:,8) = [f(2:Nr-1,:,8);f(Nr-1:Nr,:,8)]; % comp vertical
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % restoring the cells that have crossed the static boundary (bounce back)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Bounce-back (no-slip)
    G = f; 
    G(vec1) = f(vec3);
    G(vec2) = f(vec4);
    G(vec3) = f(vec1);
    G(vec4) = f(vec2);
    G(vec5) = f(vec7);
    G(vec6) = f(vec8);
    G(vec7) = f(vec5);
    G(vec8) = f(vec6);
    
    % Specular reflexion (free-Slip)
    G(vec12) = f(vec12);
    G(vec22) = f(vec42);
    G(vec32) = f(vec32);
    G(vec42) = f(vec22);
    G(vec52) = f(vec82);
    G(vec62) = f(vec72);
    G(vec72) = f(vec62);
    G(vec82) = f(vec52);
    f = G;
    
    % Block 5.2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % recalculating rho and u
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rho=sum(f,3);
    
    % Determining the velocities according to Eq.() (see slides)    
    ux = (C_x(1).*f(:,:,1)+C_x(2).*f(:,:,2)+C_x(3).*f(:,:,3)+C_x(4).*f(:,:,4)+C_x(5).*f(:,:,5)+C_x(6).*f(:,:,6)+C_x(7).*f(:,:,7)+C_x(8).*f(:,:,8))./rho ;
    uy = (C_y(1).*f(:,:,1)+C_y(2).*f(:,:,2)+C_y(3).*f(:,:,3)+C_y(4).*f(:,:,4)+C_y(5).*f(:,:,5)+C_y(6).*f(:,:,6)+C_y(7).*f(:,:,7)+C_y(8).*f(:,:,8))./rho ;
    
    % keeping the results
    rho_inlet_t(ta) = mean(rho(2:20,60)-rho_l, 'all');
    rho_outlet_t(ta) = mean(rho(2:20,251)-rho_l, 'all');
    ux_outlet_t(ta) = mean(ux(2:20,251), 'all');

    
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
    
    % Zhou et al. (2009) Axisymmetric Condition    
    dUy_dx=gradient(uy,1,2);
    dUy_dy=(1./(rho.*tau.*cs2)).*((f(:,:,2)-feq(:,:,2))-(f(:,:,4)-feq(:,:,4))+(f(:,:,5)-feq(:,:,5))+(f(:,:,6)-feq(:,:,6)) -(f(:,:,7)-feq(:,:,7)) -(f(:,:,8)-feq(:,:,8)))-(uy./(2*R));
    dUx_dy=(-1./(rho.*tau.*cs2)).*((f(:,:,5)-feq(:,:,5))-(f(:,:,6)-feq(:,:,6))+(f(:,:,7)-feq(:,:,7))-(f(:,:,8)-feq(:,:,8)))-dUy_dx;
    
    Fx = rho.*((-ux.*uy./R)+(nu.*dUx_dy./R));
    Fy = rho.*((-uy.^2./R)+(nu.*dUy_dy./R)-(nu.*uy./(R.^2)));
    
    
    for i = 1:9
        h_1(:,:,i) = -rho.*(uy./(9.*R)) ;
        h_2(:,:,i) = (1/(3*sum(C_y.^2)*cs2))*(C_x(i).*Fx + C_y(i).*Fy);
    end
   
    % Block 5.4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Collision (relaxation) step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    f = (1-omega)*f + omega*feq...
                     - sigRight.*(feq - f_t)...       % Right ABS layer sourse term
                     - sigLeft.*(feq - f_t)...        % Left ABS layer sourse term
                     - sigTop.*(feq - f_t)...         % Top ABS layer sourse term
                     - sigInlet.*(feq - feq_inlet)... % Inlet ABS layer and velocity BC sourse term
                     + h_1 + h_2;                     % Reis et al. (2007) Axisymmetric sourse terms

    

%     % Ploting the results in real time
%     figure(3)
%     surf(rho-1), view(2), shading flat, axis equal, caxis([-.01 .01]);
%     xlim([1 Mc]);
%     ylim([1 Nr]);
%     % figure(2)
%     % plot(ux(150,:))
%     % figure(3)
%     % surf(sqrt(ux.^2+uy.^2)), view(2), shading flat, axis equal, caxis([-.01 .01])
%     %     pause(.0001)    

end %  End main time Evolution Loop
delete(prog); % close waitbar
end_time = now;
%% Finishing messages
fprintf("End time: "+datestr(end_time, 'yyyy-mm-dd HH:MM:SS')+"\n");
fprintf("Simulation Elapsed Time: "+datestr(end_time-start_time, 'HH:MM:SS')+"\n");
sim_speed = ta/((end_time-start_time)*24*60*60);
fprintf("Simulation speed: %0.2f frames per second \n",sim_speed);

% save raw data
p_outlet_t = rho_outlet_t.*cs2;
writematrix([time' p_outlet_t' ux_outlet_t'],"pressure_list3_Lucas_Schroeder_3_Zhou_"+datestr(end_time, 'yyyy-mm-dd-HHMMSS')+".csv");
%%
figure(4)
surf(ux), view(2), shading flat, axis equal, colorbar; % caxis([-0.001 0.004])
hold on;
plot(x_cylynder,y_cylinder,"k");
plot(x_buffer,y_buffer,"--k");
plot(x_inlet,y_inlet,"--k");
hold off;
colormap jet;
title("Last timestep u_x countour plot");
xlim([1 Mc]);
ylim([1 Nr]);


%% Plot Pressure history
figure(5)
plot(time,p_outlet_t);hold on;
plot(time,rho_inlet_t.*cs2);
plot(time,ux_outlet_t);hold off;
    title("Outlet Pressure");
    xlabel("Time, timesteps");
    legend("Outlet Pressure","Inlet Pressure","Outlet u_x");
    set(gca,"FontSize",20);
    ylim([-0.06 0.06]);
    xlim([0 total_time]);

%% FFT and Impedance
p_outlet_t2 = p_outlet_t;
ux_outlet_t2 = ux_outlet_t;
% p_outlet_t2(ta+1:(ta*10)) = 0;  % zero-padding
% ux_outlet_t2(ta+1:(ta*10)) = 0;

T = length(p_outlet_t2);  % Sampling duration, timesteps
Lf = length(p_outlet_t2); % Number of sampling points
Df = 1/T;                 % Frequency increment
Fs = T/Lf;                % Sampling frequency (ends up being the max possible (Fs=1) because 
ka = (linspace(0,1,Lf))*r_a*2*pi/cs;          % Helmholtz number

P_f = fft(p_outlet_t2)/Lf;
Ux = fft(ux_outlet_t2)/Lf;

Z = P_f./Ux;
Z_0 = rho_l/cs;
Ref = abs((Z./Z_0 - 1)./(Z./Z_0 + 1));

figure(7)
plot(ka,Ref,"x");
    title("Reflection Coeficient R");
    xlabel("Frequency in Helmholtz number, ka");
    ylabel("|R|");
    xlim([0 1.8]);
    ylim([0.4 1]);
    set(gca,"FontSize",20);

figure(8)
subplot(2,1,1)
plot(ka,real(Z./Z_0));
    title("Normalized Impedance");
    ylabel("Resistence");
    xlim([0 1.8]);
%     ylim([0.4 1]);
    set(gca,"FontSize",20);
subplot(2,1,2)
plot(ka,imag(Z./Z_0));
    xlabel("Frequency in Helmholtz number, ka");
    ylabel("Reactance");
    xlim([0 1.8]);
%     ylim([0.4 1]);
    set(gca,"FontSize",20);

% saves the processed data
R3_Zhou = [ka.' (Z./Z_0).' Ref.'];
save("R3_Zhou.mat","R3_Zhou");

%% Verify Walls
domain = zeros(Nr,Mc,N_c);
walls = ones(Nr,Mc,N_c);
walls(vec1)=0;
walls(vec2)=0;
walls(vec3)=0;
walls(vec4)=0;
walls(vec5)=0;
walls(vec6)=0;
walls(vec7)=0;
walls(vec8)=0;

walls(vec12)=0;
walls(vec22)=0;
walls(vec32)=0;
walls(vec42)=0;
walls(vec52)=0;
walls(vec62)=0;
walls(vec72)=0;
walls(vec82)=0;

for i=1:9
    walls(:,:,1) = and(walls(:,:,1), walls(:,:,i));
end

figure(11)
surf(walls(:,:,1)), hold on, view(2), shading flat, axis equal, caxis([0 1]), colorbar;
    title("Rigid Walls Implemeted With Crossing");
    colormap gray;
    xlim([1 Mc]);
    ylim([1 Nr]);


