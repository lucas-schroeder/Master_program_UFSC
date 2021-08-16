%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This is a simple implementation of the LBGK model.
%%% By Andrey R. da Silva, August 2010
%%% Edited by Lucas Schroeder, March 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
close all

fprintf("Strongly Coupled String\n")
start_time = now;
fprintf("Start Time: "+datestr(start_time, 'yyyy-mm-dd HH:MM:SS')+"\n");

seed = 0;
real_time_plot = 1;


% Block 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lattice size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nr = 202;                    % Number of lines   (cells in the y direction)
Mc = 202;                    % Number of columns (cells in the x direction)


% Block 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Physical parameters (macro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_p = 340;                    % Sound velocity on the fluid [m/s]
rho_p = 1.2;                  % physical density [kg/m^3]
nu_STP_p = 1.5e-5;            % Air STP Kinematic Viscosity [m^2/s]
Lx = 0.302;                   % Maximun dimenssion in the x direction [m]
Ly = 0.152;                   % Maximun dimenssion on th y direction  [m]


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
    load("seed_bouzidi.mat"); % load seed results from a previous simulation
    feq = zeros(Nr,Mc,N_c);
    fprintf("Seeding from previous results (seed_bouzidi.mat)\n")
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
f_wall = zeros(Nr,Mc,9);                       % Constructing the target dist. functions
sigRight  = zeros(Nr,Mc,9);
sigLeft   = zeros(Nr,Mc,9);
sigTop    = zeros(Nr,Mc,9);
sigBot    = zeros(Nr,Mc,9);
fs = zeros(Nr,Mc,9);

% INLET - targets for the ABC condition
ux_in = 0; % FIX ME to include flow
uy_in = 0;
rho_in = 1;

% source distribution functions (uses the Maxweall formula to fill f_t
% based on target values fo rho and u
uxsq = ux_in.^2; 
uysq = uy_in.^2; 
usq = uxsq + uysq;
f_wall(:,:,1)= rho_in*w1 .*(1 +f1*ux_in +f2.*uxsq -f3*usq);
f_wall(:,:,2)= rho_in*w1 .*(1 +f1*uy_in +f2*uysq -f3*usq);
f_wall(:,:,3)= rho_in*w1 .*(1 -f1*ux_in +f2*uxsq -f3*usq);
f_wall(:,:,4)= rho_in*w1 .*(1 -f1*uy_in +f2*uysq -f3*usq);
f_wall(:,:,5)= rho_in*w2 .*(1 +f1*(+ux_in+uy_in) +f2*(+ux_in+uy_in).^2 -f3.*usq);
f_wall(:,:,6)= rho_in*w2 .*(1 +f1*(-ux_in+uy_in) +f2*(-ux_in+uy_in).^2 -f3.*usq);
f_wall(:,:,7)= rho_in*w2 .*(1 +f1*(-ux_in-uy_in) +f2*(-ux_in-uy_in).^2 -f3.*usq);
f_wall(:,:,8)= rho_in*w2 .*(1 +f1*(+ux_in-uy_in) +f2*(+ux_in-uy_in).^2 -f3.*usq);
f_wall(:,:,9)= rho_in*w0 .*(1 - f3*usq);

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
        sigTop(Nr-D:Nr,et,i) = sigma_m*(linspace(0,D,D+1)/(D)).^2;
    end
end

% BOTTOM - left side of the domain
for i=1:9
    for et = 1:Mc    
        sigBot(1:D+1,et,i) = fliplr(sigma_m*(linspace(0,D,D+1)/(D)).^2);
    end
end




% Plot ABS Layers
figure(1)
layer=3;
arr = cat(3, sigRight(:,:,layer), sigLeft(:,:,layer), sigTop(:,:,layer), sigBot(:,:,layer));
ABS_layer = max(arr,[],3);
pcolor(ABS_layer), view(2), shading flat, colorbar; hold on;
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
total_time = 50000;% # of flow passes over the channel
transient = 0;
time = 0:1:(total_time-1);

% Rope initail conditions
w = [0 2/25:2/25:2 (-52/75+8/3):-2/75:0]; % Initial position
wm = w;
F = wm*0;
uw = zeros(1,Mc);                         % Rope initial vertical velocity
P = .333;                                 % Rope tension
rho_rope = 160;                           % Rope density
w_time = zeros(1,total_time);             % Store history position for a point in the rope
d_x = 1;
d_t = 1;
w_history = zeros(total_time, length(w));

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

    wp = ((P*gradient(gradient(w,d_x),d_x) + F)*d_t.^2)/rho_rope - wm + 2*w;
    wm = w; 
    w = wp;    
    w(1) = 0;                  % Mantain fixed ends
    w(end) = 0;
    uw(50:150) = (w-wm)./d_t;  % Compute vertical velocity
    w_time(ta) = w(51);        % Storing central position
    
    % Compute new "walls" for new position of the rope
    for iDir = 1:8
        [dir(iDir).dist, dir(iDir).col, dir(iDir).line] = crossing(Nr,Mc,iDir, 50:1:150, w+100);
    end
    
    G = f;
    Dirs = [3 4 1 2 7 8 5 6 9];
    for iDir = 1:8
        for iDist = 1:length(dir(iDir).dist)
            q = dir(iDir).dist(iDist);
            if q < 1/2
                G(dir(iDir).line(iDist),dir(iDir).col(iDist),Dirs(iDir)) =  q*(2*q+1)*f(dir(iDir).line(iDist),dir(iDir).col(iDist),iDir)...
                    +(1-2*q)*(1+2*q)*f(dir(iDir).line(iDist)-C_y(iDir),dir(iDir).col(iDist)-C_x(iDir),iDir)...
                    -q*(1-2*q)*f(dir(iDir).line(iDist)-2*C_y(iDir),dir(iDir).col(iDist)-2*C_x(iDir),iDir)...
                    + 6*W(iDir)*(C_y(iDir)*uw(dir(iDir).col(iDist)));
            elseif q >= 1/2
                G(dir(iDir).line(iDist),dir(iDir).col(iDist),Dirs(iDir)) = ...
                    (1/(q*(2*q+1)))*f(dir(iDir).line(iDist),dir(iDir).col(iDist),iDir)...
                    +((2*q-1)/q)*f(dir(iDir).line(iDist),dir(iDir).col(iDist),Dirs(iDir))...
                    +((1-2*q)/(2*q+1))*f(dir(iDir).line(iDist)+C_y(iDir),dir(iDir).col(iDist)+C_x(iDir),Dirs(iDir))...
                    +((6*W(iDir))/(q*(2*q+1)))*(C_y(iDir)*uw(dir(iDir).col(iDist)));
            end
        end
    end
    
    f = G;  
    
    
    % Block 5.2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % recalculating rho, u and F
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rho=sum(f,3);
    
    % Determining the velocities according to Eq.() (see slides)    
    ux = (C_x(1).*f(:,:,1)+C_x(2).*f(:,:,2)+C_x(3).*f(:,:,3)+C_x(4).*f(:,:,4)+C_x(5).*f(:,:,5)+C_x(6).*f(:,:,6)+C_x(7).*f(:,:,7)+C_x(8).*f(:,:,8))./rho ;
    uy = (C_y(1).*f(:,:,1)+C_y(2).*f(:,:,2)+C_y(3).*f(:,:,3)+C_y(4).*f(:,:,4)+C_y(5).*f(:,:,5)+C_y(6).*f(:,:,6)+C_y(7).*f(:,:,7)+C_y(8).*f(:,:,8))./rho ;
    
    % Compute the fluid forces upon the string
    for iForce = 1:101
        F(iForce) = (rho(iForce+49,ceil(w(iForce))+100) - rho(iForce+49,floor(w(iForce))+100))*cs2;
    end
    
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
                     - sigRight.*(feq - f_wall)...       % Right ABS layer sourse term
                     - sigLeft.*(feq - f_wall)...        % Left ABS layer sourse term
                     - sigTop.*(feq - f_wall)...         % Right ABS layer sourse term
                     - sigBot.*(feq - f_wall);           % Left ABS layer sourse term
    

    % Ploting the results in real time
    if and(mod(ta,100)==0, real_time_plot)
        figure(3)
        pcolor(sqrt(usq)), view(2), shading flat, axis equal, colorbar; hold on;
        plot(50:150,w+100,'k','LineWidth',2); hold off;
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
writematrix([time' w_time'],"list4_Lucas_Schroeder_strongly_coupled_"+datestr(end_time, 'yyyy-mm-dd-HHMMSS')+".csv");

save("last_timestep_strongly_coupled.mat","f");

%% Plot velocity profile
figure(4)
pcolor(sqrt(usq)), view(2), shading flat, axis equal,  colorbar; % caxis([-0.001 0.004])
hold on;
plot([Mc-D Mc-D],[D Nr-D],"--k", "LineWidth",2);
plot([D D],[D Nr-D],"--k", "LineWidth",2);
plot([D Mc-D],[Nr-D Nr-D],"--k", "LineWidth",2);
plot([D Mc-D],[D D],"--k", "LineWidth",2);
plot(50:150,w+100,'k','LineWidth',2); hold off;
hold off;
colormap jet;
title("Last Timestep Velocity Magnitude Countour plot");
xlim([1 Mc]);
ylim([1 Nr]);
xlabel("x position");
ylabel("y position");
set(gca,"FontSize",20);

%% Plot rope displacement history
figure(6)
plot(w_time)
    title("Rope Center Point Position History");
    xlabel("x position");
    ylabel("Displacement");
    set(gca,"FontSize",20);
    xlim([1 total_time]);


