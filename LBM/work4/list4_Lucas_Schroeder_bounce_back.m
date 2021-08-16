%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This is a simple implementation of the LBGK model.
%%% By Andrey R. da Silva, August 2010
%%% Edited by Lucas Schroeder, March 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
close all

fprintf("Bounce-bakc walls\n")
start_time = now;
fprintf("Start Time: "+datestr(start_time, 'yyyy-mm-dd HH:MM:SS')+"\n");

seed = 1;
real_time_plot = 0;
flow_pass = 4;

fprintf("Number of Flow Passes: %0.0f \n", flow_pass)

% Block 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lattice size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nr = 152;                    % Number of lines   (cells in the y direction)
Mc = 302;                    % Number of columns (cells in the x direction)


% Block 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Physical parameters (macro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_p = 340;                    % Sound velocity on the fluid [m/s]
rho_p = 1.2;                  % physical density [kg/m^3]
nu_STP_p = 1.5e-5;            % Air STP Kinematic Viscosity [m^2/s]
Lx = 0.302;                   % Maximun dimenssion in the x direction [m]
Ly = 0.152;                   % Maximun dimenssion on th y direction  [m]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% definition of the solid static boundaries by vetors using "crossing"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bottom wall
x_bot = [0 Mc-1] + 0.5;
y_bot = [1 1] + 0.5;
[vec1,vec2,vec3,vec4,vec5,vec6,vec7,vec8] = crossing3(Nr,Mc,x_bot,y_bot);

% Top wall
x_top = [0 Mc-1] + 0.5;
y_top = [Nr-1 Nr-1] + 0.5;
[vec12,vec22,vec32,vec42,vec52,vec62,vec72,vec82] = crossing3(Nr,Mc,x_top,y_top);

% Cylinder contour
r_a = 15;         % cylinder radius
C = [80 75]; % center of circle
th = 0:2*pi/1000:2*pi;
x_cyn = C(1) - r_a*cos(th);
y_cyn = C(2) - r_a*sin(th);
[vec13,vec23,vec33,vec43,vec53,vec63,vec73,vec83] = crossing3(Nr,Mc,x_cyn,y_cyn);

% Block 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lattice parameters (micro - lattice unities)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cs = 1/sqrt(3);                % lattice speed of sound
cs2 = cs^2;                    % Squared speed of sound cl^2
zeta = c_p/cs;                 % Conversion factor between physiscal and lattice units
Dx = Lx/(Mc-2);                % Lattice space (pitch) [m/lattice_length]
Dt = Dx/zeta;                  % time step [s/lattice_time]

omega = 1.88;                  % Relaxation frequency
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
if seed
    load("seed_bounce_back.mat"); % load seed results from a previous simulation
    feq = zeros(Nr,Mc,N_c);
    fprintf("Seeding from previous results (seed_bounce_back.mat)\n")
else
    f = zeros(Nr,Mc,N_c);
end


% Filling the initial distribution function (at t=0) with initial values
if sum(f.^2,"all")==0
    f(:,:,:) = rho_l/9;   
    ux = zeros(Nr, Mc);
    uy = zeros(Nr, Mc);
else
    rho = sum(f,3);
    ux = (C_x(1).*f(:,:,1)+C_x(2).*f(:,:,2)+C_x(3).*f(:,:,3)+C_x(4).*f(:,:,4)+C_x(5).*f(:,:,5)+C_x(6).*f(:,:,6)+C_x(7).*f(:,:,7)+C_x(8).*f(:,:,8))./rho ;
    uy = (C_y(1).*f(:,:,1)+C_y(2).*f(:,:,2)+C_y(3).*f(:,:,3)+C_y(4).*f(:,:,4)+C_y(5).*f(:,:,5)+C_y(6).*f(:,:,6)+C_y(7).*f(:,:,7)+C_y(8).*f(:,:,8))./rho ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% ABC parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = 30;                                     % thickness of the ABC layer
sigma_m = 0.3;                              % constant - empirical
f_in = zeros(Nr,Mc,9);                       % Constructing the target dist. functions
sigRight  = zeros(Nr,Mc,9);
sigLeft   = zeros(Nr,Mc,9);
sigTop    = zeros(Nr,Mc,9);
sigInlet  = zeros(Nr,Mc,9);
fs = zeros(Nr,Mc,9);

% INLET - targets for the ABC condition
ux_in = 0.1*cs; % FIX ME to include flow
uy_in = 0;
rho_in = 1;

% source distribution functions (uses the Maxweall formula to fill f_t
% based on target values fo rho and u
uxsq = ux_in.^2; 
uysq = uy_in.^2; 
usq = uxsq + uysq;
f_in(:,:,1)= rho_in*w1 .*(1 +f1*ux_in +f2.*uxsq -f3*usq);
f_in(:,:,2)= rho_in*w1 .*(1 +f1*uy_in +f2*uysq -f3*usq);
f_in(:,:,3)= rho_in*w1 .*(1 -f1*ux_in +f2*uxsq -f3*usq);
f_in(:,:,4)= rho_in*w1 .*(1 -f1*uy_in +f2*uysq -f3*usq);
f_in(:,:,5)= rho_in*w2 .*(1 +f1*(+ux_in+uy_in) +f2*(+ux_in+uy_in).^2 -f3.*usq);
f_in(:,:,6)= rho_in*w2 .*(1 +f1*(-ux_in+uy_in) +f2*(-ux_in+uy_in).^2 -f3.*usq);
f_in(:,:,7)= rho_in*w2 .*(1 +f1*(-ux_in-uy_in) +f2*(-ux_in-uy_in).^2 -f3.*usq);
f_in(:,:,8)= rho_in*w2 .*(1 +f1*(+ux_in-uy_in) +f2*(+ux_in-uy_in).^2 -f3.*usq);
f_in(:,:,9)= rho_in*w0 .*(1 - f3*usq);

% INLET - right side of the domain
for i=1:9   
    for et = 1:Nr    
    sigRight(et,Mc-D:Mc,i) = sigma_m*(linspace(0,D,D+1)/(D)).^2;
    end
end

% OUTLET - targets for the ABC condition
f_out = f_in;

% OUTLET - left side of the domain
for i=1:9
    for et = 1:Nr    
    sigLeft(et,1:D+1,i) = fliplr(sigma_m*(linspace(0,D,D+1)/(D)).^2);
    end
end

% Plot ABS Layers
figure(1)
layer=3;
ABS = max(sigRight(:,:,layer), sigLeft(:,:,layer));
pcolor(ABS), view(2), shading flat, colorbar; hold on;
plot(x_cyn,y_cyn,"k", "LineWidth",2);
plot(x_top,y_top,"k", "LineWidth",2);
plot(x_bot,y_bot,"k", "LineWidth",2);
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
total_time = ceil(flow_pass*Mc/ux_in); % # of flow passes over the channel
transient = ceil(0*Mc/ux_in);
time = 0:1:(total_time-1);

% Generate measurements vector
rho_outlet_t = zeros(1,total_time);
rho_inlet_t = zeros(1,total_time);
ux_outlet_t= zeros(1,total_time);
ux_inlet_t= zeros(1,total_time);
rho_center_t= zeros(1,total_time);
ux_center_t= zeros(1,total_time);

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

    G = f;
    % Specular reflexion (free-Slip)
    G(vec1) = f(vec1);
    G(vec2) = f(vec4);
    G(vec3) = f(vec3);
    G(vec4) = f(vec2);
    G(vec5) = f(vec8);
    G(vec6) = f(vec7);
    G(vec7) = f(vec6);
    G(vec8) = f(vec5);
    
    G(vec12) = f(vec12);
    G(vec22) = f(vec42);
    G(vec32) = f(vec32);
    G(vec42) = f(vec22);
    G(vec52) = f(vec82);
    G(vec62) = f(vec72);
    G(vec72) = f(vec62);
    G(vec82) = f(vec52);
    
    % no-slip
    G(vec13) = f(vec33);
    G(vec23) = f(vec43);
    G(vec33) = f(vec13);
    G(vec43) = f(vec23);
    G(vec53) = f(vec73);
    G(vec63) = f(vec83);
    G(vec73) = f(vec53);
    G(vec83) = f(vec63);
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
    rho_inlet_t(ta) = mean(rho(2:Nr-2,D+2)-rho_l, 'all');
    ux_inlet_t(ta)  = mean(ux(2:Nr-2,D+2), 'all');
    rho_outlet_t(ta)= mean(rho(2:Nr-2,Mc-D-2)-rho_l, 'all');
    ux_outlet_t(ta) = mean(ux(2:Nr-2,Mc-D-2), 'all');
    rho_center_t(ta)= mean(rho(2:Nr-2,75+100)-rho_l, 'all');
    ux_center_t(ta) = mean(ux(2:Nr-2,75+100), 'all');

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
    
    
   
    % Block 5.4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Collision (relaxation) step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    f = (1-omega)*f + omega*feq...
                     - sigRight.*(feq - f_out)...       % Right ABS layer sourse term
                     - sigLeft.*(feq - f_in);          % Left ABS layer sourse term
    
 
    % Ploting the results in real time
    if and(mod(ta,200)==0, real_time_plot)
        figure(3)
        pcolor(sqrt(usq)), view(2), shading flat, axis equal; hold on;
        plot(x_cyn,y_cyn,"k", "LineWidth",2);
        plot(x_top,y_top,"k", "LineWidth",2);
        plot(x_bot,y_bot,"k", "LineWidth",2);
        plot([Mc-D Mc-D],[1 Nr],"--k", "LineWidth",2);
        plot([D D],[1 Nr],"--k", "LineWidth",2); hold off;
        colormap jet;
        xlim([1 Mc]);
        ylim([1 Nr]);
        % figure(2)
        % plot(ux(150,:))
        % figure(3)
        % surf(sqrt(ux.^2+uy.^2)), view(2), shading flat, axis equal, caxis([-.01 .01])
        % pause(.0001)
    end


end %  End main time Evolution Loop
delete(prog); % close waitbar
end_time = now;
%% Finishing messages
fprintf("End time: "+datestr(end_time, 'yyyy-mm-dd HH:MM:SS')+"\n");
fprintf("Simulation Elapsed Time: "+datestr(end_time-start_time, 'HH:MM:SS')+"\n");
sim_speed = ta/((end_time-start_time)*24*60*60);
fprintf("Simulation speed: %0.2f frames per second \n",sim_speed);
fprintf("Simulated time: %0.0f timesteps \n",ta);

% save raw data
p_outlet_t = rho_outlet_t.*cs2;
p_inlet_t = rho_inlet_t.*cs2;
p_center_t = rho_center_t.*cs2;
writematrix([time' p_inlet_t' ux_inlet_t' p_outlet_t' ux_outlet_t' p_center_t' ux_center_t'],"pressure_list4_Lucas_Schroeder_bounce_back_"+datestr(end_time, 'yyyy-mm-dd-HHMMSS')+".csv");

save("last_timestep.mat","f");

%% Plot velocity profile
figure(4)
pcolor(sqrt(usq)), view(2), shading flat, axis equal,  colorbar; % caxis([-0.001 0.004])
hold on;
plot(x_cyn,y_cyn,"k", "LineWidth",2);
plot(x_top,y_top,"k", "LineWidth",2);
plot(x_bot,y_bot,"k", "LineWidth",2);
plot([Mc-D Mc-D],[1 Nr],"--k", "LineWidth",2);
plot([D D],[1 Nr],"--k", "LineWidth",2);
hold off;
colormap jet;
title("Last Timestep Velocity Magnitude Countour plot");
xlim([1 Mc]);
ylim([1 Nr]);
xlabel("x position");
ylabel("y position");
set(gca,"FontSize",20);

%% Plot Pressure history
figure(5)
subplot(2,1,1)
    plot(time,p_outlet_t);hold on;
    plot(time,p_center_t);
    plot(time,p_inlet_t);hold off;
    title("Pressure and Velocity History");    
    legend("p_{outlet}","p_{center}","p_{inlet}");
    ylabel("Pressure, latice units");
    set(gca,"FontSize",20);
    ylim([-2 2]*1e-3);
    xlim([0 total_time]);
subplot(2,1,2)
    plot(time,ux_outlet_t);hold on;
    plot(time,ux_center_t);
    plot(time,ux_inlet_t);hold off;
    legend("Outlet u_x","Center u_x","Inlet u_x");
    xlabel("Time, timesteps");
    ylabel("Velocity, latice units");
    xlim([0 total_time]);
    set(gca,"FontSize",20);

%% Strouhal
T = length(rho_center_t);  % Sampling duration, timesteps
Lf = length(rho_center_t); % Number of sampling points
Df = 1/T;                  % Frequency increment
Fs = T/Lf;                 % Sampling frequency (ends up being the max possible (Fs=1) because we measure every timestep
freq = linspace(0,Fs,Lf);

rho_freq = fft(rho_center_t)/Lf; % Scale for number of samples

figure(6)
loglog(freq, abs(rho_freq),"-x");
    xlim([0 0.5]);

[peak, index] = max(abs(rho_freq(2:end)));
char_freq = freq(index+1)

St = char_freq*2*r_a/mean(ux_inlet_t)

