%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This is a simple implementation of the LBGK model.
%%% By Andrey R. da Silva, August 2010
%%% Edited by Lucas Schroeder, February 2021
%%% The code does not take into account any specific boundary condition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, clc, close all

fprintf("Case numbers:\n")
fprintf("1 - closed-closed\n")
fprintf("2 - closed-open\n")
fprintf("3 - closed-open with PML\n")
fprintf("4 - closed-open, with open space\n")
n = input('Enter a case number: ');

switch n
    case 1
        disp('You have chosen Case 1 - closed-closed')
    case 2
        disp('You have chosen Case 2 - closed-open')
    case 3
        disp('You have chosen Case 3 - closed-open with PML')
    case 4
        disp('You have chosen Case 4 - closed-open, with open space')
    otherwise
        disp('You have entered a invalid case')
        return
end

start_time = now;
fprintf("Start Time: "+datestr(start_time, 'yyyy-mm-dd HH:MM:SS')+"\n");


% Block 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lattice size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nr = 62;                    % Number of lines   (cells in the y direction)
% Mc = 252;                   % Number of columns (cells in the x direction)
switch n
    case {1,2,3}
        Nr = 62;                    % Number of lines   (cells in the y direction)
        Mc = 252;                   % Number of columns (cells in the x direction)
    case 4 % extra room for wave
        Nr = 186;%62;                    % Number of lines   (cells in the y direction)
        Mc = 350;%252;                   % Number of columns (cells in the x direction)
end


% Block 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Physical parameters (macro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_p = 340;                    % Sound velocity on the fluid [m/s]
rho_p = 1.2;                  % physical density [kg/m^3]
nu_STP_p = 1.5e-5;            % Air STP Kinematic Viscosity [m^2/s]
Lx = .3;                      % Maximun dimenssion in the x direction [m]
Ly = 0.072;                   % Maximun dimenssion on th y direction  [m]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% definition of the solid static boundaries by vetors using "crossing"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Tube walls
switch n
    case 1 % closed-closed, four walls
        x_keypoints = [251.5 1.5  1.5 251.5 251.5];
        y_keypoints = [  1.5 1.5 61.5  61.5   1.5];
    case {2,3} % closed-open, three walls
        x_keypoints = [251.5 1.5  1.5 251.5];
        y_keypoints = [  1.5 1.5 61.5  61.5];
    case 4 % closed-open, lift up
        x_keypoints = [251.5 1.5  1.5 251.5];
        y_keypoints = [  1.5 1.5 61.5  61.5] + 62;
end

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ABC parameters (non-reflecting boundary)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if or(n==3, n==4) % case 2 and 3 - implement absorbing layer
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


    for i=1:9                                   % reconstructing f_t with feq. Notice t    
    for et = 1:Nr    
        sigr(et,Mc-D:Mc,i) = sigma_m*(linspace(0,D,D+1)/(D)).^2;
    end
    end

end
% 
%  for i=1:9                                   % reconstructing f_t with feq. Notice t    
%     for et = 1:Nr    
%     sigl(et,1:D+1,i) = fliplr(sigma_m*(linspace(0,D,D+1)/(D)).^2);
%     end
%     end
% 


%Block 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin the iteractive process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

total_time = ceil(20*Mc/cs);
time = 0:1:(total_time-1);
pressure_history = zeros(1,total_time);
ux_history = zeros(1,total_time);

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
    
   
    % Application of the initial disturbance at t = 0;
    if ta == 1
    switch n
        case {1,2,3}
        rho(2:61,2) = 1 + .1;%*sin(2*pi*freq*(ta-1));
        case 4
        rho((2:61)+62,2) = 1 + .1;%*sin(2*pi*freq*(ta-1));
    end    
    end

    % imposicao do grdiente de pressao para barreira aberta

    %     rho(2:Nr-1,2:3) = 1+.01;%*sin(2*pi*freq*(ta-1));
    %     rho(2:Nr-1,Nr-2:Nr-1) = 1-.01;%*sin(2*pi*freq*(ta-1));
   
          
    rt0= w0*rho;
    rt1= w1*rho;
    rt2= w2*rho;

    % Determining the velocities according to Eq.() (see slides)    
    ux = (C_x(1).*f(:,:,1)+C_x(2).*f(:,:,2)+C_x(3).*f(:,:,3)+C_x(4).*f(:,:,4)+C_x(5).*f(:,:,5)+C_x(6).*f(:,:,6)+C_x(7).*f(:,:,7)+C_x(8).*f(:,:,8))./rho ;
    uy = (C_y(1).*f(:,:,1)+C_y(2).*f(:,:,2)+C_y(3).*f(:,:,3)+C_y(4).*f(:,:,4)+C_y(5).*f(:,:,5)+C_y(6).*f(:,:,6)+C_y(7).*f(:,:,7)+C_y(8).*f(:,:,8))./rho ;
    
    % keeping the results
    switch n
        case {1, 2}
            pressure_history(ta) = mean(rho(4:58,248)-1, 'all')*cs2; % isothermal model [lattice uinits]
            ux_history(ta) = mean(ux(4:58,248), 'all');
        case 3 % before absorbing layer
            pressure_history(ta) = mean(rho(4:58,218)-1, 'all')*cs2; % isothermal model [lattice uinits]
            ux_history(ta) = mean(ux(4:58,218), 'all');
        case 4 % end of tube
            pressure_history(ta) = mean(rho((4:58)+62,248)-1, 'all')*cs2; % isothermal model [lattice uinits]
            ux_history(ta) = mean(ux((4:58)+62,248), 'all');
    end
    
    % Block 5.3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determining the relaxation functions for each direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
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
   
    switch n
        case {1,2}
            f = (1-omega)*f + omega*feq;% - sigr.*(feq-f_t);% - sigl.*(feq-f_t);;
        case {3,4}%3 % absorbing layer change the source term in LBGK
            f = (1-omega)*f + omega*feq - sigr.*(feq - f_t);
    end
    
 
%     % Ploting the results in real time  \
%     figure(1)
%     surf(rho-1), view(2), shading flat, axis equal, caxis([-.01 .01])
%     hold on;
%     plot(x_keypoints,y_keypoints,"k");
%     hold off;
%     % surf(rho-1), view(2), shading flat, axis equal, caxis([-.001 .001])
%     % colormap(gray)
%     % grid on
%     % figure(2)
%     % plot(ux(150,:))
%     % figure(3)
%     % surf(sqrt(ux.^2+uy.^2)), view(2), shading flat, axis equal, caxis([-.01 .01])
%     % pause(.0001)    

end %  End main time Evolution Loop
delete(prog);
end_time = now;
fprintf("End time: "+datestr(end_time, 'yyyy-mm-dd HH:MM:SS')+"\n");
fprintf("Simulation Elapsed Time: "+datestr(end_time-start_time, 'HH:MM:SS')+"\n");

%% Plot results
nfig = 1;
figure(nfig)
plot(time,ux_history);
hold on;
plot(time,pressure_history);hold off;
title("Acoustic particle velocity, component x");
xlabel("Time, timesteps");
ylabel("Ux, latticeVelocity");
legend("particle velocity","acoustic pressure");

nfig = nfig + 1;
figure(nfig)
surf(rho-1), view(2), shading flat, axis equal, caxis([-.05 .05]), colorbar;
hold on;
plot(x_keypoints,y_keypoints,"k");
% plot(248,(4:58)+62,"kx");
hold off;
title("Last timestep density countour plot");
grid on
colormap jet;
xlim([1 Mc]);
ylim([1 Nr]);

%% Calculate impedance

P_f = fft(pressure_history);
U_f = fft(ux_history);
L = length(pressure_history);
Z0 = rho_l*cs;
Z = P_f./U_f;
Df = 1/Dt;
% time_fft = Dt*(0:(L - 1));
frq = Df*(0:(L - 1)/2)/L;
kh = 2*pi*frq*Ly./c_p;

nfig = nfig + 1;
figure(nfig)
subplot(2,1,1)
    plot(kh, real(Z(1:length(kh))),"o");
    title("Impedance at tube's termination");
    ylabel("Normalized Resistence");
    grid on;
    xlim([0 2]);
%     ylim([0.55 0.8]);
subplot(2,1,2)
    plot(kh, imag(Z(1:length(kh))),"o");
    ylabel("Normalized Reactance");
    xlabel("Helmholtz number");
    grid on;
    xlim([0 2]);
%     ylim([-0.05 0.1]);
    
%% cuting the results
% time2 = time(350:750);
% ux_history2 = ux_history(350:750);
% pressure_history2 = pressure_history(350:750);
% 
% P = fft(pressure_history2);
% U = fft(ux_history2);
% Z0 = rho_l*cs;
% Z = P./U;
% Df = 1/Dt;
% time2 = Dt*(0:(length(pressure_history2)-1));
% frequency = Df*(0:(length(pressure_history2)-1)/2)/length(pressure_history2);
% helmholtz = 2*pi*frequency*Ly./c_p;
% 
% figure(7)
% subplot(2,1,1)
%     plot(helmholtz, real(Z(1:length(helmholtz)))/Z0);
%     ylabel("Resistence");
%     xlim([0 15]);
%     ylim([-40 10]);
% subplot(2,1,2)
%     plot(helmholtz, imag(Z(1:length(helmholtz)))/Z0);
%     ylabel("Reactance");
%     xlim([0 15]);
%     ylim([-40 10]);

%% another way to calculate impedance
% time_p = time.*Dt;
% ux_history_p = ux_history*zeta;
% pressure_history_p = pressure_history*rho_p*c_p^2;
% 
% 
% Fs = 1/Dt;                    % Sampling frequency
% L = length(pressure_history); % Length of signal
% fq = Fs*(0:(L/2))/L;
% 
% [Z_, frq] = tfestimate(ux_history, pressure_history,hamming(L),0,L,Fs);
% helmholtz = 2*pi*frq*Ly./c_p;
% Z0=rho_l*cs;
% 
% % Z_ = Z_(frq <= (c_p/(2*Ly))); % only frequencies below cut-off
% % helmholtz = helmholtz(frq <= (c_p/(2*Ly)));
% figure(5)
% subplot(2,1,1)
%     plot(helmholtz, real(Z_)/Z0,"o");
%     ylabel("Resistence");
%     xlim([0 15]);
%     ylim([-5 20]);
% subplot(2,1,2)
%     plot(helmholtz, imag(Z_)/Z0,"o");
%     ylabel("Reactance");
%     xlim([0 15]);
%     ylim([-60 5]);
%% yet another way to calculate impedance


% aux1 = fft(pressure_history)/L;
% aux2 = aux1(1:L/2+1);
% P_f = 2*aux2;
% 
% aux3 = fft(ux_history)/L;
% aux4 = aux3(1:L/2+1);
% Ux_f = 2*aux4;
% 
% figure(3)
% plot(fq, abs(P_f));
% 
% Z = P_f./Ux_f;
% figure(4)
% plot(fq, abs(Z));

% pn = (rho(150,150:300) - 1)*(1/sqrt(3))^2;
% xn = linspace(1,151,151);
% A = .001/(cs^2);
% 
% [p,x] = cylin_wave(freq,visc,cs,A,linspace(1,149,149),0);
% 
% plot(x,p), hold on
% plot(xn,pn)
