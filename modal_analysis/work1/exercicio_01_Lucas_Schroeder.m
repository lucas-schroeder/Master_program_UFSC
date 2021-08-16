clc; clear;
%% System parameters

m = 1;      % [kg]
k = 1.58e6; % [N/m]
eta = 0.05; % [-]

%% 1 - Excitation

F = 1; % [N], excitation amplitude
range_start = 100;  % [Hz]
range_finish = 500; % [Hz]
df = 0.5;           % [Hz], frequency increment
f = range_start:df:range_finish;
w = 2*pi*f;

receptance = @(w) 1./(k - w.^2*m + 1j*k*eta); 

H = receptance(w);
H = H.*(1 + 1e-2*randn(1,length(f)));

nfig = 1;
figure(nfig);
semilogx(f, abs(H));
title("Receptance Magnitude");
xlabel("");
ylabel("");

nfig = nfig + 1;
figure(nfig);
semilogx(f, real(H));
title("Receptance - Real Part");
xlabel("Frequency, Hz");
ylabel("Re(H), m/N");

nfig = nfig + 1;
figure(nfig);
semilogx(f, imag(H));
title("Receptance - Imaginary Part");
xlabel("Frequency, Hz");
ylabel("Im(H), m/N");

nfig = nfig + 1;
figure(nfig);
plot(real(H), imag(H));
title("Receptance - Nyquist Plot");
xlabel("Re(H), m/N");
ylabel("Im(H), m/N");

nfig = nfig + 1;
figure(nfig);
plot3(f, real(H), imag(H), "o");
title("Receptance - 3D Plot");
xlabel("Frequency, Hz");
ylabel("Re(H), m/N");
zlabel("Im(H), m/N");

%% 1.1 Find the system's damping using the Halfpower Method

[pks,central_frq] = findpeaks(abs(H),f, 'MinPeakHeight',0.1e-5,...
                                        'MinPeakDistance',200);
drop = pks*10^(-3/20); % 3dB drop in signal strength

FRQ = @(frf) interp1(abs(H),f,frf,"spline"); % interpolates and returns the frequency
bandwidth = abs(2*(central_frq - FRQ(drop)))
eta_estimated = bandwidth/central_frq

%% 1.2 Find the bandwidth of the Halfpower Method using only the real part

[pk1,f_a] = findpeaks(real(H),f, 'MinPeakHeight',2e-6,...
                                 'MinPeakDistance',200);
[pk2,f_b] = findpeaks(-real(H),f, 'MinPeakHeight',2e-6,...
                                 'MinPeakDistance',200);
central_frq2 = (f_a+f_b)/2;
bandwidth2 = abs(f_b-f_a)
eta_estimated2 = abs(f_b-f_a)/central_frq2

%% 1.3 Find the halp power point using the Nyquist diagram

[pk1,f_a] = findpeaks(real(H),f, 'MinPeakHeight',2e-6,...
                                 'MinPeakDistance',200);
[pk2,f_b] = findpeaks(-real(H),f, 'MinPeakHeight',2e-6,...
                                 'MinPeakDistance',200);
[pk3,central_frq3] = findpeaks(-imag(H),f, 'MinPeakHeight',2e-6,...
                                          'MinPeakDistance',200);
                                      
bandwidth3 = abs(f_b-f_a)
eta_estimated3 = bandwidth3/central_frq3

%% 1.4 Find the system's damping using the Imaginary part of the Receptance


%% 1.5 Find the dissipated energy in each cycle at ressonance

%%
FRF = @(frq) interp1(f,H,frq,"spline"); % Interpolation of the FRF