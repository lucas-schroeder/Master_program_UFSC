% Lucas Schroeder, Feb. 25, 2021
clear all; close all; clc

%% Parameter
N = 3; % # of degrees of fredom
[m1, m2, m3] = deal(1., .95, 1.05);
[k1, k2, k3, k4, k5, k6] = deal(1e3, 1e3, 1e3, 1e3, 1e3, 1e3);


M = diag([m1 m2 m3]);
K = [k1+k4+k6   -k4      -k6   ;
       -k4    k2+k4+k5   -k5   ;
       -k6      -k5    k3+k5+k6];
   

%% 1.1 Proportional Structural Damping

beta = 0.05;
gamma = 0.0;
D = beta.*K + gamma.*M;
[Psi, W] = eig((K+1j.*D),M);
wr = diag(real(W).^0.5);
fn = wr./(2*pi);                   % Natural Frequencies, Hz
eta_r = beta + gamma./(wr.^2);     % Modal damping
csi_r = eta_r./2;
Lambda2 = real(W).*(1 + 1j.*eta_r);   % Equal to W, just to check
m_r = diag(transpose(Psi)*M*Psi);     % Modal masses
m_r_aux = repmat(transpose(m_r),N,1); % prepare for element-wise multiplication
Phi = Psi.*(m_r_aux.^-0.5);           % Nomalize eigenvectors


%% mode shapes
% figure('Name', 'Mode Shapes','NumberTitle','off');
% set(gca,'OuterPosition',[6.3543e-01  -1.3311e-01   2.7454e-01   1.1932e+00]);
% subplot(1,3,1);
%     plot([0; Phi(:,1)],0:N,'-o');
%     title('First Mode','Interpreter','latex','FontSize',12);
%     ylabel('DoF','Interpreter','latex','FontSize',16); 
%     ylim([0.99 4]);
%     xlabel('Amplitude','Interpreter','latex','FontSize',16);
%     xlim([-1 1]);
%     set(gca,'FontSize',14);
% subplot(1,3,2);
%     plot([0; Phi(:,2)],0:N,'-o');
%     title('Second Mode','Interpreter','latex','FontSize',12);
%     xlabel('Amplitude','Interpreter','latex','FontSize',16);  
%     ylim([0.99 4]);
%     xlim([-1 1]);
% %     ylim([-0.15 0.15]);
%     set(gca,'FontSize',14);
% subplot(1,3,3);
%     plot([0; Phi(:,3)],0:N,'-o');
%     title('Third Mode','Interpreter','latex','FontSize',12);
%     xlabel('Amplitude','Interpreter','latex','FontSize',16);
%     ylim([0.99 4]);
%     xlim([-1 1]);
% %     ylim([-0.15 0.15]);
%     set(gca,'FontSize',14);
% print -deps fig_modos

% Mode shapes in complex plane
figure('Name', 'Mode Shapes - Complex','NumberTitle','off');
set(gca,'OuterPosition',[6.3543e-01  -1.3311e-01   2.7454e-01   1.1932e+00]);
subplot(1,3,1);
    compass(Phi(:,1));
    title('First Mode','Interpreter','latex','FontSize',12);
    ylabel('Im(Phi) - First mode','Interpreter','latex','FontSize',16); 
    xlabel('Re(Phi) - First mode','Interpreter','latex','FontSize',16);
%     ylim([-1 1]);
%     xlim([-1 1]);
    set(gca,'FontSize',14);
subplot(1,3,2);
    compass(Phi(:,2));
    title('Second Mode','Interpreter','latex','FontSize',12);
%     ylabel('Im(Phi) - Second mode','Interpreter','latex','FontSize',16); 
    xlabel('Re(Phi) - Second mode','Interpreter','latex','FontSize',16);
%     ylim([-1 1]);
%     xlim([-1 1]);
    set(gca,'FontSize',14);
subplot(1,3,3);
    compass(Phi(:,3));
    title('Third Mode','Interpreter','latex','FontSize',12);
%     ylabel('Im(Phi) - Third mode','Interpreter','latex','FontSize',16); 
    xlabel('Re(Phi) - Second mode','Interpreter','latex','FontSize',16);
%     ylim([-1 1]);
%     xlim([-1 1]);
    set(gca,'FontSize',14);
print -deps fig_modos

%% 1.2 Receptance calculation
Fs = 20;              % Upper limit frequency
df = 0.01;            % Frequency increment FIX ME
L = Fs/df;            % Length of frequencies vector
f = df*(0:L-1);       % Frequencies vector
w = 2*pi.*f;
dt = 1/Fs;            % Time increment
t = dt*(0:L-1);       % Time vector

k = 1;       % DoF of incedent force

H = zeros(max(size(w)),3); % Receptância do GL j em relação ao GL k
A_original=zeros(N,N,N);   % Modal constants 3d matrix
for j=1:3                % for each response Dof
    for p=1:max(size(w)) % through all the frequencies
        for n=1:N        % sum for each mode n
            A_original(j,k,n) = Phi(j,n)*Phi(k,n);
            H(p,j) = H(p,j) +             A_original(j,k,n)/...
                              (wr(n).^2  - w(p).^2 + 2j*csi_r(n)*wr(n)*w(p));
            
        end
    end
end

%% 1.2 Plot receptance

nfig=1;
figure(nfig)
for k=1:N % Response of k-th DoF
subplot(2,N,k)
    semilogy(f, abs(H(:,k)),"x");
    title(sprintf("Receptance H_{%d1}",k));
    xlim([0 Fs])
    % ylim([0.05 100])
    ylabel(sprintf("|H_{%d1}|",k))
    xlabel("Frequency, Hz")
    grid on;
    set(gca,'FontSize',18);
subplot(2,N,k+N)
    plot(f, wrapToPi(phase(H(:,k))));
    xlim([0 Fs])
    ylim([-4 4])
    ylabel(sprintf("Phase of H_{%d1}",k))
    xlabel("Frequency, Hz")
    grid on;
    set(gca,'FontSize',18);
end

k=1;
nfig=nfig+1;
figure(nfig)
semilogy(f, abs(H(:,k)),"x");
    title(sprintf("Receptance H_{%d1}",k));
    xlim([9.8 10.4])
    ylim([2.4 2.6].*1e-3)
    ylabel(sprintf("|H_{%d1}|",k))
    xlabel("Frequency, Hz")
    grid on;
    set(gca,'FontSize',18);

%% 1.3 Plot Nyquist
nfig=nfig+1;
figure(nfig)
for k=1:N % Response of k-th DoF
subplot(1,N,k)
    plot(real(H(:,k)), imag(H(:,k)), "o");
    title(sprintf("Nyquist plot of H_{%d1}",k));
%     xlim([0 Fs])
%     ylim([0.05 100])
    ylabel(sprintf("Im(H_{%d1})",k))
    xlabel(sprintf("Re(H_{%d1})",k));
    grid on;
    set(gca,'FontSize',18);
end

nfig=nfig+1;
figure(nfig)
for k=1:N % Response of k-th DoF
subplot(2,N,k)
    plot(f, real(H(:,k)),"x");
    title(sprintf("Real part H_{%d1}",k));
    ylabel(sprintf("Re(H_{%d1})",k));
    xlabel("Frequency, Hz");
    grid on;
    set(gca,'FontSize',18);
subplot(2,N,k+N)
    plot(f, imag(H(:,k)),"x");
    title(sprintf("Imaginary part H_{%d1}",k));
    ylabel(sprintf("Im(H_{%d1})",k));
    xlabel("Frequency, Hz");
    grid on;
    set(gca,'FontSize',18);
end

% nfig=nfig+1;
% figure(nfig)
% subplot(2,1,1)
%     plot(f, real(H(:,k)), 'g');
%     xlim([0 Fs])
%     ylabel("Re(H_{11}), m/N");
%     grid on;
% subplot(2,1,2)
%     plot(f, imag(H(:,k)), 'b');
%     xlim([0 Fs])
%     % ylim([-40 40])
%     grid on;
%     ylabel("Im(H_{11}), m/N");
%     xlabel("Frequency, Hz")
%     
% nfig=nfig+1;
% figure(nfig);
% plot(real(H(:,k)), imag(H(:,k)), "o");
% title("Receptance - Nyquist");
% xlabel("Re(H_{11}), m/N");
% ylabel("Im(H_{11}), m/N");
%     
% nfig=nfig+1;
% figure(nfig);
% plot3(f, real(H(:,k)), imag(H(:,k)), "o");
% title("Receptance - 3D Plot");
% xlabel("Frequency, Hz");
% ylabel("Re(H_{11}), m/N");
% zlabel("Im(H_{11}), m/N");
%% 1.3 Find ressonant frequencies

for j=1:size(H,2) % DoF of response (H_jk)
[peak_at_ressonance,res_frq] = findpeaks(abs(H(:,j)),f, 'MinPeakHeight',1.1e-3,...
                                                        'MinPeakDistance',0.1);
for r=1:length(res_frq)
    % FIX ME! - change peak bandwidth seach
    around_f = and(f>=(res_frq(r)*0.9),f<=(res_frq(r)*1.1));   % only near ressonance
    around_pk = abs(H(:,j))>=peak_at_ressonance(r)/2;          % amp grater than half power (1.414)
    around = and(around_f,around_pk');
%     figure(r+10)
%     plot(f(around),abs(real(H(around,j))))
    [amp_at_ressonance,aux1] = findpeaks(abs(real(H(around,j))),f(around), 'MinPeakHeight',1.1e-3,...
                                                                           'MinPeakDistance',0.015);

    diameter = sum(amp_at_ressonance);
    Ar_j1(j,r)=diameter.*(2*pi*res_frq(r)).^2.*eta_r(r); % rAjk, modal contant for input k=1
    
    % We know the diameter, but dont konw the signal of rPhi_j.rPhi_k
    if imag(H((f==res_frq(r)),j))>0
        Ar_j1(j,r) = - Ar_j1(j,r);
    end
end
end

Ar_j1
real(A_original(:,1,:))

%% 1.3 Reconstruct Receptance
csi_rec = csi_r; % "measured" modal damping for each mode
H_rec = zeros(max(size(w)),size(Ar_j1,1)); % Receptância do GL j em relação ao GL k
N_rec = size(Ar_j1,2); % number of modes
N_dof = size(Ar_j1,1); % number of degrees of fredom
for j=1:N_dof            % for each response Dof
    for p=1:max(size(w)) % through all the frequencies
        for r=1:N_rec    % sum for each mode n
            H_rec(p,j) = H_rec(p,j) +          Ar_j1(j,r)/...
                            ((2*pi*res_frq(r)).^2  - w(p).^2 + 2j*csi_rec(r)*(2*pi*res_frq(r))*w(p));
        end
    end
end

% Plot reconstructed Receptance 
nfig=nfig+1;
figure(nfig)
for k=1:N_dof % Response of k-th DoF
subplot(2,N_dof,k)
    semilogy(f, abs(H_rec(:,k)),"x"); hold on;
    semilogy(f, abs(H(:,k))); hold off;
    title(sprintf("Receptance H_{%d1}",k));
    xlim([0 Fs])
    % ylim([0.05 100])
    legend("Reconstructed", "Original",'FontSize',12);
    set(gca,'FontSize',12);
    ylabel(sprintf("|H_{%d1}|",k))
    xlabel("Frequency, Hz")
    grid on;
    set(gca,'FontSize',18);
subplot(2,N_dof,k+N_dof)
    plot(f, wrapToPi(phase(H_rec(:,k)))); hold on;
    plot(f, wrapToPi(phase(H(:,k)))); hold off;
    xlim([0 Fs])
    ylim([-4 4])
    ylabel(sprintf("Phase of H_{%d1}",k))
    xlabel("Frequency, Hz")
    grid on;
    set(gca,'FontSize',18);
end
%% 1.6 Inverse
nfig=nfig+1;
figure(nfig)
for j=1:3 % Response of k-th DoF
subplot(3,1,j)
    plot(f,imag(1./H(:,j)));
    if j==1; title(sprintf("Inverse of receptance H_{j1}, for j=1,2,3")); end
    xlim([0 Fs])
    ylabel(sprintf("Im(1/H_{%d1})",j))    
    grid on;
    set(gca,'FontSize',18);
end
xlabel("Frequency, Hz")
