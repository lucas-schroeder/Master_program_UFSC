% Lucas Schroeder, Mar. 8, 2021
clear all; close all; clc

%% Parameter
N = 12;     % # of degrees of fredom
L = 1;      % beam length, m
dx = L/N;   % mesh resolution
x = dx/2:dx:(L-dx/2);
m = 1.2;    % beam mass, kg
eta = 0.05; % structural damping

%% 1 Eigenvector and eigenvalues

Phi = zeros(N);
for n=1:N     % each mode
    Phi(:,n) = sqrt(2/m)*sin(n*pi.*x./L);
end

fn = [67.70 270.60 608.86 1082.4];
wn = 2*pi*fn;


%% 2 Plot mode shapes
nfig=1;
figure(nfig)
for j=1:4
subplot(4,1,j)
    plot(0:N+1,[0; Phi(:,j); 0],'-o');
    grid on;
%     title('First Mode','Interpreter','latex','FontSize',12);
    xlabel('DoF','Interpreter','latex','FontSize',16);
    ylabel('Mode Amplitude','Interpreter','latex','FontSize',16);    
    xlim([0 13]);
    ylim([-1.5 1.5]);
    set(gca,'FontSize',14);    
end



%% 3 Receptance calculation
Fs = 1600;            % Upper limit frequency
df = 1;               % Frequency increment FIX ME
L = Fs/df;            % Length of frequencies vector
f = df*(0:L-1);       % Frequencies vector
w = 2*pi.*f;
dt = 1/Fs;            % Time increment
t = dt*(0:L-1);       % Time vector

k = 1;       % DoF of incedent force

H = zeros(max(size(w)),N); % Receptância do GL j em relação ao GL k
A_original=zeros(N,N,N);   % Modal constants 3D matrix
for j=1:N                % for each response Dof
    for p=1:max(size(w)) % through all the frequencies
        for n=1:4        % sum for each mode n
            A_original(j,k,n) = Phi(j,n)*Phi(k,n);
            H(p,j) = H(p,j) +             A_original(j,k,n)/...
                              (wn(n).^2  - w(p).^2 + 1j*eta*wn(n).^2);            
        end
    end
end

% H = H + (0.2.*H).*rand(size(H)); % add noise

%% 4 Plot receptance

nfig=1;
figure(nfig)
for k=1:N % Response of k-th DoF
subplot(2,N,k)
    semilogy(f, abs(H(:,k)));
    title(sprintf("Receptance H_{%d1}",k));
    xlim([0 Fs])
    % ylim([0.05 100])
    ylabel(sprintf("|H_{%d1}|",k))
    xlabel("Frequency, Hz")
    grid on;
    set(gca,'FontSize',7);
subplot(2,N,k+N)
    plot(f, wrapToPi(phase(H(:,k))));
    xlim([0 Fs])
    ylim([-4 4])
    ylabel(sprintf("Phase of H_{%d1}",k))
    xlabel("Frequency, Hz")
    grid on;
    set(gca,'FontSize',7);
end

%% 5 Nyquist Plot
nfig=nfig+1;
figure(nfig)
index = f<=400;
for k=1:N % Response of k-th DoF
subplot(6,2,k)
    plot(real(H(index,k)), imag(H(index,k)));
    title(sprintf("Nyquist plot of H_{%d1}",k));
%     xlim([-2 2].*1e-6)
%     ylim([-4 0].*1e-6)
%     ylabel(sprintf("Im(H_{%d1})",k))
%     xlabel(sprintf("Re(H_{%d1})",k));
    grid on;
    set(gca,'FontSize',8);
end

%% 6 Real and Imaginary Plot
nfig=nfig+1;
figure(nfig)
for k=1:N % Response of k-th DoF
subplot(2,N,k)
    plot(f, real(H(:,k)));
    title(sprintf("Real part H_{%d1}",k));
    ylabel(sprintf("Re(H_{%d1})",k));
    xlabel("Frequency, Hz");
    grid on;
    set(gca,'FontSize',8);
subplot(2,N,k+N)
    plot(f, imag(H(:,k)));
    title(sprintf("Imaginary part H_{%d1}",k));
    ylabel(sprintf("Im(H_{%d1})",k));
    xlabel("Frequency, Hz");
    ylim([-9e-6 5e-6])
    grid on;
    set(gca,'FontSize',8);
end


%% 7 Find Ressonant Frequencies, Modal Constants and Modal Damping

for j=1:size(H,2) % DoF of response (H_jk)
[peak_at_ressonance,res_frq] = findpeaks(abs(H(:,j)),f, 'MinPeakHeight',1e-7,...
                                                        'MinPeakDistance',100);
                                                    
eta_guess=0.1;
for r=1:length(res_frq) % mode
    bw = eta_guess*res_frq(r);
    around_f = and(f>=(res_frq(r)-bw/2),f<=(res_frq(r)+bw/2));   % only near ressonance
    around_pk = abs(H(:,j))>=peak_at_ressonance(r)/2;            % amp grater than half power (1.414)
    around = and(around_f,around_pk');
%     figure(r+10)
%     plot(f(around),abs(real(H(around,j))))
    [amp_at_ressonance,halfPowerFreq] = findpeaks(abs(real(H(around,j))),f(around), 'MinPeakHeight',1e-8,...
                                                                                    'MinPeakDistance',2);
                                                   
	HalfPowerBandwith = abs(halfPowerFreq(2)-halfPowerFreq(1));
    eta_rec_aux(j,r) = HalfPowerBandwith/res_frq(r);
    diameter(j,r) = sum(amp_at_ressonance);
    
    % We know the diameter, but dont konw the signal of rPhi_j.rPhi_k
    if imag(H((f==res_frq(r)),j))>0
        diameter(j,r) = - diameter(j,r);
    end
end
end
eta_rec = mean(eta_rec_aux,1); % mean of each column, each mode has a damping parameter for all DoF
Ar_j1 = diameter.*(2*pi*res_frq).^2.*eta_rec;
% Ar_j1
% real(A_original(:,1,:))

%% 8 Reconstruct Receptance
Fs = 1600;            % Upper limit frequency
df = 1;               % Frequency increment FIX ME
L = Fs/df;            % Length of frequencies vector
f = df*(0:L-1);       % Frequencies vector
w = 2*pi.*f;
dt = 1/Fs;            % Time increment
t = dt*(0:L-1);       % Time vector

H_rec = zeros(max(size(w)),size(Ar_j1,1)); % Receptância do GL j em relação ao GL k
N_rec = size(Ar_j1,2); % number of modes
N_dof = size(Ar_j1,1); % number of degrees of fredom
for j=1:N_dof            % for each response Dof
    for p=1:max(size(w)) % through all the frequencies
        for r=1:N_rec    % sum for each mode n
            H_rec(p,j) = H_rec(p,j) +          Ar_j1(j,r)/...
                            ((2*pi*res_frq(r)).^2 - w(p).^2 + 1j*eta_rec(r)*(2*pi*res_frq(r)).^2);
        end
    end
end

% Plot reconstructed Receptance 
nfig=nfig+1;
figure(nfig)
N_dof = 12;
for k=1:N_dof % Response of k-th DoF
subplot(2,N_dof,k)
    semilogy(f, abs(H(:,k))); hold on;
    semilogy(f, abs(H_rec(:,k))); hold off;
    title(sprintf("Receptance H_{%d1}",k),'FontSize',8);
    xlim([0 Fs])
    if k==1
    legend("Original","Reconstructed",'FontSize',10);
    ylabel(sprintf("|H_{1}|"))
    end
    grid on;
subplot(2,N_dof,k+N_dof)
    plot(f, wrapToPi(phase(H(:,k)))); hold on;
    plot(f, wrapToPi(phase(H_rec(:,k)))); hold off;
    if k==1
%     legend("Reconstructed", "Original",'FontSize',8);
    ylabel(sprintf("Phase of H_{j1}"))
    end
    xlim([0 Fs])
    ylim([-4 4])
    xlabel("Frequency, Hz")
    grid on;
end
%% 9 Reconstruct Eigenvectors from modal constants

% Ar_j1(j,r) -> % rAjk, modal contant for input k=1, output j and mode r
N_rec = size(Ar_j1,2); % number of modes
N_dof = size(Ar_j1,1); % number of degrees of fredom j
Phi2 = zeros(N_dof, N_rec);

for n=1:N_rec
Phi2(1,n) = sqrt(abs(Ar_j1(1,n))); % we dont know the signal of 1_Phi_1
Phi2(:,n) = (1/Phi2(1,n)).*Ar_j1(:,n);
[Phi2(:,n) Phi(:,n)]               % check eigenvector
end



%% 10 Plot reconstructed mode shapes
nfig=nfig+1;
figure(nfig)
for j=1:4
subplot(4,1,j)
    plot(0:N+1,[0; Phi(:,j); 0]); hold on;
    plot(0:N+1,[0; Phi2(:,j); 0],'x'); hold off;
    grid on;
    title(sprintf('Mode %0.0f',j),'Interpreter','latex','FontSize',16);
    xlabel('DoF','Interpreter','latex','FontSize',14);
    ylabel('Mode Amplitude','Interpreter','latex','FontSize',14);    
    xlim([0 13]);
    ylim([-1.5 1.5]);
%     set(gca,'FontSize',14);    
end
legend("Original", "Reconstructed",'FontSize',12);
