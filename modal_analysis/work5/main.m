%%
% By Lucas Schroeder - 02 Apr 2021

clear
clc
close all

nfig=1;

%% 1.1 Import data - Hammer

opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [77, 1101];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Freq", "Real", "Imag"];
opts.VariableTypes = ["double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data (original Y axis unit	g/N)
dof = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33];
for j=1:length(dof)
    j = uint32(dof(j));
    file_name = "FRF 3-"+int2str(j)+".txt";
    aux = readtable("Dados martelo\" + file_name, opts);
    f = aux.Freq;
    FRF(:,j) = aux.Real + 1j.*aux.Imag;
end
FRF = FRF./9.81;            % Convert g/N to (m/s²)/N
FRF = FRF./(-(2*pi*f).^2);  % Convert from Accelerance ((m/s²)/N) to Receptance (m/N)
clear opts

fn = [
    17.212
    68.840
    118.878
    155.391
    226.395
    265.219
    235.230
    362.805
    414.118
    567.539
    602.292
    643.435
    707.200
    740.567
    820.439
    880.725
    943.641
    ];
res_frq = round(fn).';
res_w = 2*pi*res_frq;

%% 1.2 Find peaks, Nyquis diameter, and modal dampings

for j=1:length(dof)    % for each dof of response (H_jk)
    for r=1:length(fn) % for each ressonant freq (mode)
        arround = 3;
        interval = and(f>round(fn(r))-arround, f<round(fn(r))+arround);
        interval_start = find(interval,1)-1;

        [pks(r,dof(j)), locs_index] = max(abs(imag(FRF(interval,dof(j))))); %#ok<*SAGROW>
        locs(r,dof(j)) = locs_index + interval_start; % only works because dF = 1 Hz

        [pks_real(r,dof(j)), f1_index] = max(real(FRF(interval,dof(j))));   
        [pks_real_(r,dof(j)), f2_index] = max(-real(FRF(interval,dof(j))));
        f1 = f1_index + interval_start;
        f2 = f2_index + interval_start;
        
        HalfPowerBandwith = abs(f2 - f1);
        diameter(j, r) = abs(pks_real(r,dof(j))) + abs(pks_real_(r,dof(j)));
        eta_rec_aux(j,r) = HalfPowerBandwith/locs(r,dof(j));
               
        % We know the diameter, but dont konw the signal of rPhi_j.rPhi_k
        if imag(FRF(locs(r,dof(j)),j))<=0
            diameter(j,r) = - diameter(j,r);
        end
        
    end
end

j = 1; % dof of response (H_jk)
r = 2; % ressonant freq (mode)

nfig=nfig+1;
figure(nfig)
plot(f(locs(r,dof(j))), abs(pks(r,dof(j))),"o", "LineWidth",2); hold on;
yline(pks_real(r,dof(j)));
yline(-pks_real_(r,dof(j)));
plot(f, imag(FRF(:,dof(j))));
plot(f, real(FRF(:,dof(j))));
plot(f, abs(FRF(:,dof(j)))); hold off;
legend("Max Magnitude", "Diameter limits", "Diameter limits", "Imag", "Real", "Abs")
set(gca,"FontSize",18)

%% 1.3 Calculate Modal constants
eta_rec = mean(eta_rec_aux,1); % mean of each column, each mode has a damping parameter for all DoF;
Ar_j2 = diameter.*(2*pi*res_frq).^2.*eta_rec;

%% 1.4 Reconstruct Eigenvectors from modal constants
N = 15; % number points in longitudinal direction

% Ar_j2(j,r) -> % rAjk, modal contant for input k=2, output j and mode r
N_mode = size(Ar_j2,2); % number of modes
N_dof = size(Ar_j2,1); % number of degrees of fredom j
Phi = zeros(N_dof, N_mode);

% Phi(j,r) -> % rPhi_j, mode shape for dof j and mode r
for n=1:N_mode
    Phi(2,n) = sqrt(abs(Ar_j2(2,n)));
    Phi(:,n) = (1/Phi(2,n)).*Ar_j2(:,n);
end

% Plot
figure(2)
for mode=1:4
    subplot(4,1,mode)
    plot(1:17,[0; Phi(1:15,mode); 0],'-x');
    xlim([1 17]);
    ylabel("\psi");
    title("Mode " + num2str(mode));
    grid on;
    set(gca,"FontSize",16);
end
xlabel("DoF");

%% 2.1 Plot real and imaginary parts of the inverse Receptance
mode = 2;
k = 3;

arround = 8;
interval = and(f>locs(mode,dof(k))-arround, f<locs(mode,dof(k))+arround);

 
subplot(2,1,1)
    plot(f(interval), real(1./FRF(interval,dof(k)))); hold on;
%     plot(locs(mode,dof(k)), real(1./FRF(f==locs(mode,dof(k)),dof(k))),"o", "LineWidth",2); hold off;
    yline(0);
    title("Mode "+num2str(mode)+" Response DoF "+num2str(k));
    ylabel("Re(H^{-1})");
    set(gca,"FontSize",16);
    grid on
subplot(2,1,2)
    plot(f(interval), imag(1./FRF(interval,dof(k)))); hold on;
%     plot(locs(mode,dof(k)), imag(1./FRF(f==locs(mode,dof(k)),dof(k))),"o", "LineWidth",2)
    yline(0);
    ylabel("Im(H^{-1})");
    xlabel("Frequncy, Hz");
    set(gca,"FontSize",16);
    grid on


%% 3.1 Reconstruct Receptance and ajust with Residual Stifness of higher modes
Fs = 1024;            % Upper limit frequency
df = 0.1;               % Frequency increment FIX ME
L = Fs/df;            % Length of frequencies vector
frq = df*(0:L-1);       % Frequencies vector
w = 2*pi.*frq;
dt = 1/Fs;            % Time increment
t = dt*(0:L-1);       % Time vector

H_rec = zeros(max(size(w)),size(Ar_j2,1)); % Receptância do GL j em relação ao GL k
% N_mode = size(Ar_j2,2); % number of modes
N_dof = size(Ar_j2,1); % number of degrees of fredom
for j=1:N_dof            % for each response Dof
    for p=1:max(size(w)) % through all the frequencies
        for r=1:4    % sum for each mode n
            H_rec(p,j) = H_rec(p,j) +          Ar_j2(j,r)/...
                            ((2*pi*res_frq(r)).^2 - w(p).^2 + 1j*eta_rec(r)*(2*pi*res_frq(r)).^2);
        end
    end
end



k = 2; % response dof
K(dof(k))=1e8;
loglog(frq, abs(H_rec(:,dof(k)) + 1/K(dof(k))),"LineWidth",1);hold on;
loglog(f, abs(FRF(:,dof(k))), "LineWidth",1); hold off
    xlim([1 1025]);
    title("Regenerated Receptance, using first four modes, DoF: "+num2str(k));
    ylabel("H(f), m/N");
    xlabel("Frequency, Hz");
    legend("Regenerated","Measured");
    grid on;
    set(gca,"FontSize",16);



%% 3.2 Plot reconstructed Receptance 
k = 2; % response dof

nfig = nfig+1;
figure(nfig)
semilogy(frq, abs(H_rec(:,dof(k))));hold on;
semilogy(f, abs(FRF(:,dof(k)))); hold off
xlim([0 1025]);


N_dof = 12;

% for k=1:N_dof % Response of k-th DoF
% subplot(2,N_dof,k)
%     semilogy(f, abs(H(:,k))); hold on;
%     semilogy(f, abs(H_rec(:,k))); hold off;
%     title(sprintf("Receptance H_{%d1}",k),'FontSize',8);
%     xlim([0 Fs])
%     if k==1
%     legend("Original","Reconstructed",'FontSize',10);
%     ylabel(sprintf("|H_{1}|"))
%     end
%     grid on;
% subplot(2,N_dof,k+N_dof)
%     plot(f, wrapToPi(phase(H(:,k)))); hold on;
%     plot(f, wrapToPi(phase(H_rec(:,k)))); hold off;
%     if k==1
% %     legend("Reconstructed", "Original",'FontSize',8);
%     ylabel(sprintf("Phase of H_{j1}"))
%     end
%     xlim([0 Fs])
%     ylim([-4 4])
%     xlabel("Frequency, Hz")
%     grid on;
% end

%% 4.1 Import data - Shaker

opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [77, 1101];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Freq", "Real", "Imag"];
opts.VariableTypes = ["double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data (original Y axis unit	g/N)
dof = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33];
for j=1:length(dof)
    j = uint32(dof(j));
    file_name = "FRF "+int2str(j)+"-2.txt";
    aux = readtable("Dados shaker ruido branco\" + file_name, opts);
    f_s = aux.Freq;
    FRF_s(:,j) = aux.Real + 1j.*aux.Imag;
end
FRF_s = FRF_s./9.81; % Convert g/N to (m/s²)/N
FRF_s = FRF_s./(-(2*pi*f_s).^2); % From Acelerance to Receptance
clear opts

%% 4.2 Fit a circle to the second mode
% Second mode is fn = 68 Hz
arround = 3;
interval = and(f_s>round(68)-arround, f_s<round(68)+arround);
interval_start = find(interval,1)-1;

for q=1:length(dof)/2
    j = dof(q);
%     nfig=nfig+1;
%     figure(nfig);
%     plot3(f_s(interval),real(FRF_s(f_s(interval),j)), imag(FRF_s(f_s(interval),j)));

    points_to_fit = [real(FRF_s(f_s(interval),j)), imag(FRF_s(f_s(interval),j))];
    Par = CircleFitByPratt(points_to_fit);
    r_a = Par(3);         % circle radius
    C = [Par(1) Par(2)];  % center of circle

    % A2_j2(j) -> % rAjk, modal contant for input k=2, output j and mode r=2
    A2_j2(j) = (2*r_a).*(2*pi*68.5).^2.*eta_rec(j);
    
    % We know the diameter, but dont konw the signal of rPhi_j.rPhi_k
    if imag(FRF_s(f_s==68,j))<=0
        A2_j2(j) = -A2_j2(j);
    end

    th = 0:2*pi/1000:2*pi;
    x_fit = C(1) - r_a*cos(th);
    y_fit = C(2) - r_a*sin(th);

    subplot(5,3,q)
    plot(real(FRF_s(f_s(interval),j)), imag(FRF_s(f_s(interval),j)), "-o"); hold on;
    plot(x_fit, y_fit);
    title("Dof "+ num2str(j));
end

%% 4.3 Reconstruct Eigenvectors from modal constants
N_dof = 16; % number points in longitudinal direction
Phi2 = zeros(N_dof,1);

% Phi2(j) -> % rPhi_j, mode shape for dof j and mode r=2
Phi2(2) = sqrt(abs(A2_j2(2)));
Phi2(:) = -(1/Phi2(2)).*A2_j2(:);

nfig=nfig+1;
figure(nfig);
plot(1:17,[Phi2(:); 0],'-x'); hold on;
plot(1:17,[0; Phi(1:15, 2); 0],'-x'); hold off;
    xlim([1 17]);
    ylabel("\psi");
    xlabel("DoF");
    title("Mode shape for Mode 2");
    legend("Shaker","Hammer");
    grid on;
    set(gca,"FontSize",16);


%% 5.1 Find analytical mode shapes and frequencies
N_anlt = 17; % Number of modes

E = 186e9;   % Modulus of elasticity, Pa
rho = 7880;  % Density, kg/m^3
L = 0.800;   % Length, m
b = 0.050;   % Width, m
h = 0.005;   % Heigth, m
m = L*b*h*rho; % Beam mass, kg
I = b*h^3/12;  % Moment of inertial, m^4


wn_anlt = (1:1:N_anlt).^2.*(pi^2).*sqrt(E*I/(m*L^3));
fn_anlt = wn_anlt./(2*pi);

x = (0:0.05:L);
Phi_anlt = zeros(length(x),N_anlt);
for n=1:N_anlt     % each mode
    Phi_anlt(:,n) = sqrt(2/m)*sin(n*pi.*x./L);
end

% figure(2)
% mode = 2;
% plot(x, Phi_anlt(:,mode),'-x'); hold on;

nfig=nfig+1;
figure(nfig);
for mode=1:4
    subplot(4,1,mode)
    plot(1:17, -Phi_anlt(:,mode)/8,'-x',"LineWidth",2); hold on;
    plot(1:17,[0; Phi(1:15,mode); 0],'-x',"LineWidth",2); hold off;
    xlim([1 17]);
    ylabel("\psi");
    title("Mode " + num2str(mode));
    grid on;
    set(gca,"FontSize",16);
end
xlabel("DoF");
legend("Regenerated","Measured");

%% 5.2 Compare Analytical and Experimental FRF
eta=0.01;

% calculate analytical receptance
Fs = 1025;            % Upper limit frequency
df = 1;               % Frequency increment FIX ME
L = Fs/df;            % Length of frequencies vector
f_anlt = df*(0:L-1);  % Frequencies vector
w_anlt = 2*pi.*f_anlt;
dt = 1/Fs;            % Time increment
t = dt*(0:L-1);       % Time vector

k = 2;       % DoF of incedent force
H = zeros(max(size(w_anlt)),N_anlt);    % Receptance of dof j to a force applied at dof k
A_original=zeros(N_anlt,N_anlt,N_anlt); % Modal constants 3D matrix
for j=1:length(x)                       % For each response Dof ...
    for p=1:max(size(w_anlt))           % through all the frequencies ...
        for n=1:N_anlt                  % sum each mode contribution.
            A_original(j,k,n) = Phi_anlt(j,n)*Phi_anlt(k,n); % nAjk, modal contant for input k, output j and mode n
            H(p,j) = H(p,j) +             A_original(j,k,n)/...
                              (wn_anlt(n).^2  - w_anlt(p).^2 + 1j*eta*wn_anlt(n).^2);            
        end
    end
end


subplot(2,1,1)
    k = 13; % Response dof
    loglog(f_anlt, abs(H(:,k)).*0.05); hold on;
    loglog(f, abs(FRF(:,k))); hold off
    ylabel("|H(f)|, m/N");
    title("Receptance of Dof " + num2str(k));
    grid on;
    set(gca,"FontSize",16);
subplot(2,1,2)
    semilogx(f_anlt, angle(H(:,k))); hold on;
    semilogx(f, angle(FRF(:,k))); hold off
    ylabel("arg(H(f)), rad");
    xlabel("Frequency, Hz");
    grid on;
    set(gca,"FontSize",16);


%% 6 Compute the Modal Assurance Criteria (MAC) of the first four flexing mode shapes
% Phi(j,r) -> % rPhi_j, mode shape for dof j and mode r
for mode_a = 1:4
    for mode_ex = 1:4
        MAC(mode_a,mode_ex) = (Phi(1:17,mode_ex).'*Phi_anlt(1:17,mode_a))^2/...
                              (Phi(1:17,mode_ex).'*Phi(1:17,mode_ex))*(Phi_anlt(1:17,mode_a).'*Phi_anlt(1:17,mode_a));
    end
end

b = bar3(MAC);
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
colorbar
title("MAC")
xlabel("Analytical Mode");
ylabel("Experimental Mode");
set(gca,"FontSize",16);

%%
for mode_a = 1:17
    for mode_ex = 1:17
        autoMAC(mode_a,mode_ex) =                (Phi(1:17,mode_ex).'*Phi(1:17,mode_a))^2/...
                               (Phi(1:17,mode_ex).'*Phi(1:17,mode_ex))*(Phi(1:17,mode_a).'*Phi(1:17,mode_a));
    end
end

%% Dynamic Vibration Absorber (DVA)
% k = 185000;
% m = 2;

k = 2*1.85e5; % N/m
m = 2;        % kg
c = 200;      % N/(m/s)
w_rec=2*pi.*frq;

H_dva = (k - m*w_rec.^2 + 1j.*w_rec*c)./((-m*w_rec.^2).*(k + 1j.*w_rec.*c));

[pk, ind] = min(abs(H_dva));
% frq(ind)

nfig=nfig+1;
figure(nfig);
loglog(frq, abs(H_dva))
    title("DVA Receptance");
    ylabel("|H_{DVA}|, m/N");
    xlabel("Frequency, Hz");
    xlim([1 1025])
    ylim([1e-7 1])
    grid on;
    set(gca,"FontSize", 22);

%%
% Compute H_jj, point receptance for DVA
j = 5; % DoF where the DVA was installed
H_jj = zeros(length(w_rec),1);
for p=1:max(size(w_rec))        % Through all the frequencies ...
    for n=1:17                  % sum each mode contribution.
        H_jj(p) = H_jj(p) +           Phi(j,n)*Phi(j,n)./... % nPj, mode shape for dof j and mode n
                          (res_w(n).^2  - w_rec(p).^2 + 1j*eta_rec(n)*res_w(n).^2);            
    end
end

% Compute H_qj, between input DVA force and output displacement
q = 2; % DoF of response
H_qj = zeros(length(w_rec),1);
for p=1:max(size(w_rec))        % Through all the frequencies ...
    for n=1:17                  % sum each mode contribution.
        H_qj(p) = H_qj(p) +           Phi(q,n)*Phi(j,n)./... % nPj, mode shape for dof j and mode n
                          (res_w(n).^2  - w_rec(p).^2 + 1j*eta_rec(n)*res_w(n).^2);            
    end
end

% Compute H_ji, between input force and DVA
i = 2; % DoF of input force
H_ji = zeros(length(w_rec),1);
for p=1:max(size(w_rec))        % Through all the frequencies ...
    for n=1:17                  % sum each mode contribution.
        H_ji(p) = H_ji(p) +           Phi(j,n)*Phi(i,n)./... % nPj, mode shape for dof j and mode n
                          (res_w(n).^2  - w_rec(p).^2 + 1j*eta_rec(n)*res_w(n).^2);            
    end
end

% Compute H_qi, between input force and output displacement
H_qi = zeros(length(w_rec),1);
for p=1:max(size(w_rec))        % Through all the frequencies ...
    for n=1:17                  % sum each mode contribution.
        H_qi(p) = H_qi(p) +           Phi(q,n)*Phi(i,n)./... % nPj, mode shape for dof j and mode n
                          (res_w(n).^2  - w_rec(p).^2 + 1j*eta_rec(n)*res_w(n).^2);            
    end
end


% H_qi_corr -> response at dof q to a force at dof i, accounting for DVA at dof j
H_qi_corr = H_qi - (H_qj.*H_ji)./(H_dva.' + H_jj);

nfig=nfig+1;
figure(nfig);
loglog(frq, abs(H_qi),"--", "LineWidth", 1.2); hold on;
loglog(frq, abs(H_qi_corr), "LineWidth", 1.4); hold off;
    title("Beam Receptance");
    legend("Original", "With ADV");
    xlim([1 1025]);
    ylabel("|H|");
    xlabel("Frequency, Hz");
    grid on;
    set(gca,"FontSize", 22);








