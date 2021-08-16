clear 
clc
close all

% Signal caracteristics
amp = 1;   % signal amplitude [-]
frq = 100; % signal frequency [Hz]
phi = 0; 
f1 = @(t) amp*sin(2*pi*frq*t-phi);


%% 1
% Sampling parameters
Fs = 400 ; % Sampling frequency [Hz]
Dt = 1/Fs; % Sampling period [s]
T = 1.0;   % Sampling duration [s]
Lf = T/Dt; % Number of sampling points
Df = Fs/Lf;% Frequency increment [Hz]
t = (0:Lf-1)*Dt; % Time vector

% Plot sampling points
nfig=1;
plot(t, f1(t),"-x")
    title("Sampled signal");
    ylabel("Amplitude");
    xlabel("Time, s");
    ylim([-1.5 1.5]);
    set(gca,"Fontsize", 20);

% FFT
F1 = fft(f1(t))/Lf;           % scale for number of sampling points
f  = linspace(0,Fs/2,Lf/2+1); % frequency vector, cut mirrored part, up to Nyquist freq
F1 = 2*F1(1:length(f));       % double the power of spectrum

% RMS
rms1 = sqrt(sum(abs(F1).^2)/2);

% Plot Spectrum
nfig=nfig+1;
figure(nfig);
subplot(2,1,1)
    stem(f, abs(F1(1:length(f))));
    title("Frequency domain signal");
    ylabel("Amplitude");
    ylim([0 2]);
    set(gca,"Fontsize", 20);
subplot(2,1,2)
    stem(f, angle(F1(1:length(f))));
    ylabel("Phase, rad");
    xlabel("Frequency, Hz");
    set(gca,"Fontsize", 20);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2
% Sampling parameters
Fs = 400 ; % Sampling frequency [Hz]
Dt = 1/Fs; % Sampling period [s]
T = 0.25;   % Sampling duration [s]
Lf = T/Dt; % Number of sampling points
Df = Fs/Lf;% Frequency increment [Hz]
t = (0:Lf-1)*Dt; % Time vector

% Plot sampling points
nfig=1;
plot(t, f1(t),"x")
    title("Sampled signal");
    ylabel("Amplitude");
    xlabel("Time, s");
    set(gca,"Fontsize", 20);

% FFT
F1 = fft(f1(t))/Lf;          % scale for number of sampling points
f  = linspace(0,Fs/2,Lf/2+1); % frequency vector, cut mirrored part, up to Nyquist freq
F1 = 2*F1(1:length(f));      % double the power of spectrum

% RMS
rms2 = sqrt(sum(abs(F1).^2)/2);

% Plot Spectrum
nfig=nfig+1;
figure(nfig);
subplot(2,1,1)
    stem(f, abs(F1(1:length(f))));
    title("Frequency domain signal");
    ylabel("Amplitude");
    ylim([0 2]);
    set(gca,"Fontsize", 20);
subplot(2,1,2)
    stem(f, angle(F1(1:length(f))));
    ylabel("Phase, rad");
    xlabel("Frequency, Hz");
    set(gca,"Fontsize", 20);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3
% Sampling parameters
Fs = 400 ; % Sampling frequency [Hz]
Dt = 1/Fs; % Sampling period [s]
T = 0.255;   % Sampling duration [s]
Lf = T/Dt; % Number of sampling points
Df = Fs/Lf;% Frequency increment [Hz]
t = (0:Lf-1)*Dt; % Time vector

% Plot sampling points
nfig=1;
plot(t, f1(t),"x")
    title("Sampled signal");
    ylabel("Amplitude");
    xlabel("Time, s");
    set(gca,"Fontsize", 20);

% FFT
F1 = fft(f1(t))/Lf;           % scale for number of sampling points
f  = linspace(0,Fs/2,Lf/2+1); % frequency vector, cut mirrored part, up to Nyquist freq
F1 = 2*F1(1:length(f));       % double the power of spectrum

% RMS
rms3 = sqrt(sum(abs(F1).^2)/2);

% Half-power Bandwidth
[pks,locs] = findpeaks(abs(F1));  % find the peak
above_peak = (f>=f(locs));        % helps find the half-power frq from de decay after the peak
half_power_frq = interp1(abs(F1(above_peak)), f(above_peak), pks/sqrt(2)); % interpolate to find the frq
bw = 2*(half_power_frq - f(locs)); % compute the bandwidth
eta3 = bw/f(locs);                 % estimate the damping

% Plot Spectrum
nfig=nfig+1;
figure(nfig);
subplot(2,1,1)
    stem(f, abs(F1));
    title("Frequency domain signal");
    ylabel("Amplitude");
    ylim([0 1]);
    set(gca,"Fontsize", 20);
subplot(2,1,2)
    stem(f, angle(F1(1:length(f))));
    ylabel("Phase, rad");
    xlabel("Frequency, Hz");
    set(gca,"Fontsize", 20);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4
% Sampling parameters
Fs = 400 ; % Sampling frequency [Hz]
Dt = 1/Fs; % Sampling period [s]
T = 0.25;   % Sampling duration [s]
Lf = T/Dt; % Number of sampling points
Df = Fs/Lf;% Frequency increment [Hz]
t = (0:Lf-1)*Dt; % Time vector

% Zero-padding
signal = [f1(t) zeros(1, 300)];
Lf = length(signal); % Number of sampling points
t = (0:Lf-1)*Dt;     % update time vector

% Plot sampling points
nfig=1;
plot(t, signal, "-x")
    title("Sampled signal with zero-padding");
    ylabel("Amplitude");
    xlabel("Time, s");
    ylim([-1.5 1.5]);
    set(gca,"Fontsize", 20);

% FFT
F1 = fft(signal)/Lf;           % scale for number of sampling points
f  = linspace(0,Fs/2,Lf/2+1);  % frequency vector, cut mirrored part, up to Nyquist freq
F1 = 2*F1(1:length(f));        % double the power of spectrum

% RMS
rms4 = sqrt(sum(abs(F1).^2)/2);
rms4_t = rms(signal);           % calculate the rms value in the time domain

% Half-power Bandwidth
[pks,locs] = findpeaks(abs(F1), 'minpeakheight', 0.1);  % find the peak
above_peak = (f>=f(locs));         % helps find the half-power frq from de decay after the peak
half_power_frq = interp1(abs(F1(above_peak)), f(above_peak), pks/sqrt(2)); % interpolate to find the frq
bw = 2*(half_power_frq - f(locs)); % compute the bandwidth
eta4 = bw/f(locs);                 % estimate the damping

% Plot Spectrum
nfig=nfig+1;
figure(nfig);
subplot(2,1,1)
    stem(f, abs(F1(1:length(f))));
    title("Frequency domain signal");
    ylabel("Amplitude");
    set(gca,"Fontsize", 20);
subplot(2,1,2)
    stem(f, angle(F1(1:length(f))));
    ylabel("Phase, rad");
    xlabel("Frequency, Hz");
    set(gca,"Fontsize", 20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5
% Signal caracteristics
amp = 1;   % signal amplitude [-]
frq = 20; % signal frequency [Hz]
f1 = @(t) amp*sin(2*pi*frq*t);

% Sampling parameters
Fs = 400 ; % Sampling frequency [Hz]
Dt = 1/Fs; % Sampling period [s]
T = 20*0.255;   % Sampling duration [s]
Lf = T/Dt; % Number of sampling points
Df = Fs/Lf;% Frequency increment [Hz]
t = (0:Lf-1)*Dt; % Time vector

% Plot sampling points
nfig=1;
plot(t, f1(t),"-x")
    title("Sampled signal");
    ylabel("Amplitude");
    xlabel("Time, s");
    xlim([0 T]);
    ylim([-1.5 1.5]);
    set(gca,"Fontsize", 20);

% FFT
F1 = fft(f1(t))/Lf;           % scale for number of sampling points
f  = linspace(0,Fs/2,Lf/2+1); % frequency vector, cut mirrored part, up to Nyquist freq
F1 = 2*F1(1:length(f));       % double the power of spectrum

% RMS
rms5 = sqrt(sum(abs(F1).^2)/2);

% Half-power Bandwidth
[pks,locs] = findpeaks(abs(F1), 'minpeakheight', 0.1);  % find the peak
above_peak = (f>=f(locs));         % helps find the half-power frq from de decay after the peak
half_power_frq = interp1(abs(F1(above_peak)), f(above_peak), pks/sqrt(2)); % interpolate to find the frq
bw = 2*(half_power_frq - f(locs)); % compute the bandwidth
eta5 = bw/f(locs);                 % estimate the damping

% Plot Spectrum
nfig=nfig+1;
figure(nfig);
subplot(2,1,1)
    stem(f, abs(F1(1:length(f))));
    title("Frequency domain signal");
    ylabel("Amplitude");
    ylim([0 1.5]);
    set(gca,"Fontsize", 20);
subplot(2,1,2)
    stem(f, angle(F1(1:length(f))));
    ylabel("Phase, rad");
    xlabel("Frequency, Hz");
    ylim([-pi pi]);
    set(gca,"Fontsize", 20);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6
% Signal caracteristics
amp = 1;   % signal amplitude [-]
frq = 100; % signal frequency [Hz]
f1 = @(t) amp*sin(2*pi*frq*t);

% Sampling parameters
Fs = 400 ; % Sampling frequency [Hz]
Dt = 1/Fs; % Sampling period [s]
T = 0.255;   % Sampling duration [s]
Lf = T/Dt; % Number of sampling points
Df = Fs/Lf;% Frequency increment [Hz]
t = (0:Lf-1)*Dt; % Time vector

% Hann Window
windowed = hann(Lf)'.*f1(t);

ECF6 = 1/rms(hann(Lf));  % Window Energy Correction Factor
% ECF = sqrt(8/3);         % Continuous Hann Window
ACF6 = 1/mean(hann(Lf)); % Window Amplitude Correction Factor

% Plot window shape
nfig=1;
figure(nfig);
plot(t,hann(Lf));
    title("Hann Window");
    ylabel("Amplitude");
    xlabel("Time, s");
    xlim([0 T]);
    set(gca,"Fontsize", 20);

% Plot sampling points
nfig=nfig+1;
figure(nfig);
plot(t, windowed,"-x")
    title("Sampled signal with Hann Window");
    ylabel("Amplitude");
    xlabel("Time, s");
    xlim([0 T]);
    ylim([-1.5 1.5]);
    set(gca,"Fontsize", 20);

% FFT
F1 = fft(windowed)/Lf;  % scale for number of sampling points
f  = linspace(0,Fs/2,Lf/2+1);   % frequency vector, cut mirrored part, up to Nyquist freq
F1 = 2*F1(1:length(f));         % double the power of spectrum
F1 = ECF6.*F1;             % scale for hann window

% RMS
rms6 = sqrt(sum(abs(F1).^2)/2);

% Half-power Bandwidth
[pks,locs] = findpeaks(abs(F1), 'minpeakheight', 0.1);  % find the peak
above_peak = (f>=f(locs));        % helps find the half-power frq from de decay after the peak
half_power_frq = interp1(abs(F1(above_peak)), f(above_peak), pks/sqrt(2)); % interpolate to find the frq
bw = 2*(half_power_frq - f(locs)); % compute the bandwidth
eta6 = bw/f(locs);                  % estimate the damping

% Plot Spectrum
nfig=nfig+1;
figure(nfig);
subplot(2,1,1)
    stem(f, abs(F1(1:length(f))));
    title("Frequency domain signal");
    ylabel("Amplitude");
    ylim([0 1]);
    set(gca,"Fontsize", 20);
subplot(2,1,2)
    stem(f, angle(F1(1:length(f))));
    ylabel("Phase, rad");
    xlabel("Frequency, Hz");
    ylim([-pi pi]);
    set(gca,"Fontsize", 20);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7
% Signal caracteristics
amp = 1;   % signal amplitude [-]
frq = 100; % signal frequency [Hz]
f1 = @(t) amp*sin(2*pi*frq*t);

% Sampling parameters
Fs = 400 ; % Sampling frequency [Hz]
Dt = 1/Fs; % Sampling period [s]
T = 1.006;   % Sampling duration [s]
Lf = T/Dt; % Number of sampling points
Df = Fs/Lf;% Frequency increment [Hz]
t = (0:Lf-1)*Dt; % Time vector

% Exponential Window
zeta = 0.03;
expwin = exp(-zeta*2*pi*frq*t);
windowed = expwin.*f1(t);

ECF7 = 1/rms(expwin);  % Window Energy Correction Factor
ACF7 = 1/mean(expwin); % Window Amplitude Correction Factor

% Plot window shape
nfig=1;
figure(nfig);
plot(t,expwin);
    title("Exponential Window");
    ylabel("Amplitude");
    xlabel("Time, s");
    xlim([0 T]);
    set(gca,"Fontsize", 20);

% Plot sampling points
nfig=nfig+1;
figure(nfig);
plot(t, windowed,"-x")
    title("Sampled signal with Exponential Window");
    ylabel("Amplitude");
    xlabel("Time, s");
    xlim([0 T]);
    set(gca,"Fontsize", 20);   

% FFT
F1 = fft(windowed)/Lf;  % scale for number of sampling points
f  = linspace(0,Fs/2,Lf/2+1);   % frequency vector, cut mirrored part, up to Nyquist freq
F1 = 2*F1(1:length(f));         % double the power of spectrum
F1 = ECF7.*F1;                  % scale for window

% RMS
rms7 = sqrt(sum(abs(F1).^2)/2);

% Half-power Bandwidth
[pks,locs] = findpeaks(abs(F1), 'minpeakheight', 0.1);  % find the peak
above_peak = (f>=f(locs));         % helps find the half-power frq from de decay after the peak
half_power_frq = interp1(abs(F1(above_peak)), f(above_peak), pks/sqrt(2)); % interpolate to find the frq
bw = 2*(half_power_frq - f(locs)); % compute the bandwidth
eta7 = bw/f(locs);                 % estimate the damping

% Plot Spectrum
nfig=nfig+1;
figure(nfig);
subplot(2,1,1)
    stem(f, abs(F1(1:length(f))));
    title("Frequency domain signal");
    ylabel("Amplitude");
    ylim([0 1]);
    set(gca,"Fontsize", 20);
subplot(2,1,2)
    stem(f, angle(F1(1:length(f))));
    ylabel("Phase, rad");
    xlabel("Frequency, Hz");
    set(gca,"Fontsize", 20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8 - Without Window
% Signal caracteristics
amp = 1;    % signal amplitude [-]
frq = 100;  % signal frequency [Hz]
eta = 0.01; % hysteretic damping
f1 = @(t) amp*sin(2*pi*frq*t).*exp(-eta*2*pi*frq*t);

% Sampling parameters
Fs = 400 ; % Sampling frequency [Hz]
Dt = 1/Fs; % Sampling period [s]
T = 1.005;   % Sampling duration [s]
Lf = T/Dt; % Number of sampling points
Df = Fs/Lf;% Frequency increment [Hz]
t = (0:Lf-1)*Dt; % Time vector

% Plot sampling points
nfig=1;
figure(nfig);
plot(t, f1(t),"-x")
    title("Sampled Signal with Hysteretic Damping");
    ylabel("Amplitude");
    xlabel("Time, s");
    xlim([0 T]);
    set(gca,"Fontsize", 20);   

% FFT
F1 = fft(f1(t))/Lf;             % scale for number of sampling points
f  = linspace(0,Fs/2,Lf/2+1);   % frequency vector, cut mirrored part, up to Nyquist freq
F1 = 2*F1(1:length(f));         % double the power of spectrum

% RMS
rms8 = sqrt(sum(abs(F1).^2)/2);

% Half-power Bandwidth
[pks,locs] = findpeaks(abs(F1), 'minpeakheight', 0.1);  % find the peak
above_peak = (f>=f(locs));         % helps find the half-power frq from de decay after the peak
half_power_frq = interp1(abs(F1(above_peak)), f(above_peak), pks/sqrt(2)); % interpolate to find the frq
bw = 2*(half_power_frq - f(locs)); % compute the bandwidth
eta8 = bw/f(locs);                 % estimate the damping

% Plot Spectrum
nfig=nfig+1;
figure(nfig);
subplot(2,1,1)
    stem(f, abs(F1(1:length(f))));
    title("Frequency domain signal");
    ylabel("Amplitude");
    ylim([0 0.2]);
    set(gca,"Fontsize", 20);
subplot(2,1,2)
    stem(f, angle(F1(1:length(f))));
    ylabel("Phase, rad");
    xlabel("Frequency, Hz");
    set(gca,"Fontsize", 20);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8 - With Window
% Signal caracteristics
amp = 1;    % signal amplitude [-]
frq = 100;  % signal frequency [Hz]
eta = 0.01; % hysteretic damping
f1 = @(t) amp*sin(2*pi*frq*t).*exp(-eta*2*pi*frq*t);

% Sampling parameters
Fs = 400 ; % Sampling frequency [Hz]
Dt = 1/Fs; % Sampling period [s]
T = 1.005;   % Sampling duration [s]
Lf = T/Dt; % Number of sampling points
Df = Fs/Lf;% Frequency increment [Hz]
t = (0:Lf-1)*Dt; % Time vector

% Exponential Window
zeta = 0.03;
expwin = exp(-zeta*2*pi*frq*t);
windowed = expwin.*f1(t);

ECF8 = 1/rms(expwin);  % Window Energy Correction Factor
ACF8 = 1/mean(expwin); % Window Amplitude Correction Factor

% Plot window shape
nfig=1;
figure(nfig);
plot(t,expwin);
    title("Exponential Window");
    ylabel("Amplitude");
    xlabel("Time, s");
    xlim([0 T]);
    set(gca,"Fontsize", 20);

% Plot sampling points
nfig=nfig+1;
figure(nfig);
plot(t, windowed,"-x")
    title("Sampled signal with Exponential Window");
    ylabel("Amplitude");
    xlabel("Time, s");
    xlim([0 T]);
    set(gca,"Fontsize", 20);   

% FFT
F1 = fft(windowed)/Lf;          % scale for number of sampling points
f  = linspace(0,Fs/2,Lf/2+1);   % frequency vector, cut mirrored part, up to Nyquist freq
F1 = 2*F1(1:length(f));         % double the power of spectrum
F1 = ECF8.*F1;                  % scale for window

% RMS
rms8b = sqrt(sum(abs(F1).^2)/2);
rms8b_t = rms(f1(t));           % calculate the rms value in the time domain

% Half-power Bandwidth
[pks,locs] = findpeaks(abs(F1), 'minpeakheight', 0.1);  % find the peak
above_peak = (f>=f(locs));         % helps find the half-power frq from de decay after the peak
half_power_frq = interp1(abs(F1(above_peak)), f(above_peak), pks/sqrt(2)); % interpolate to find the frq
bw = 2*(half_power_frq - f(locs)); % compute the bandwidth
eta8b = bw/f(locs);                 % estimate the damping

% Plot Spectrum
nfig=nfig+1;
figure(nfig);
subplot(2,1,1)
    stem(f, abs(F1(1:length(f))));
    title("Frequency domain signal");
    ylabel("Amplitude");
    ylim([0 0.3]);
    set(gca,"Fontsize", 20);
subplot(2,1,2)
    stem(f, angle(F1(1:length(f))));
    ylabel("Phase, rad");
    xlabel("Frequency, Hz");
    set(gca,"Fontsize", 20);

