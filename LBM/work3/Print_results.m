%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This code loads the results from all simulations and plots them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, clc, close all

%% Load results
load("LevineSchwinger_solucao\LevineSchwinger.mat");

mat = dir("results\*.mat");    % (struc) lists all .mat files in the folder
for ii = 1:length(mat)
    load("results\" + mat(ii).name); % uses the struc name atribute to load all .mat files
end
nfig=1;
%% 1.1 - 2D Channel
figure(nfig);
plot(LevineSchwinger(:,1), LevineSchwinger(:,3), 'LineWidth',2); hold on;
plot(R1_2DChannel(:,1), R1_2DChannel(:,3), 'LineWidth',2); hold off;
    title("2D Channel - No Flow");
    legend("Levine and Schwinger","2D Channel - No Flow")
    xlabel("Frequency in Helmholtz number, ka");
    ylabel("|R|");
    xlim([0 1.8]);
    ylim([0.4 1.1]);
    set(gca,"FontSize",16);

%% 1.2 - 2D Channel with FLOW
nfig=nfig+1;
figure(nfig);
plot(LevineSchwinger(:,1), LevineSchwinger(:,3), 'LineWidth',2); hold on;
plot(R1_2DChannel_Flow(:,1), R1_2DChannel_Flow(:,3),"x", 'LineWidth',2); hold off;
    title("2D Channel - Flow");
    legend("Levine and Schwinger","2D Channel - Flow")
    xlabel("Frequency in Helmholtz number, ka");
    ylabel("|R|");
    xlim([0 1.8]);
    ylim([0.4 1.1]);
    set(gca,"FontSize",16);

%% 2.1 - Reis Axisymmetric
nfig=nfig+1;
figure(nfig);
plot(LevineSchwinger(:,1), LevineSchwinger(:,3), 'LineWidth',2); hold on;
plot(R2_Reis(:,1), R2_Reis(:,3), "x", 'LineWidth',2); hold off;
    title("Reis et al. Axyssimetric - No Flow");
    legend("Levine and Schwinger","Reis et al.")
    xlabel("Frequency in Helmholtz number, ka");
    ylabel("|R|");
    xlim([0 1.8]);
    ylim([0.4 1.1]);
    set(gca,"FontSize",16);


%% 2.2 - Reis Axisymmetric with FLOW
nfig=nfig+1;
figure(nfig);
plot(LevineSchwinger(:,1), LevineSchwinger(:,3), 'LineWidth',2); hold on;
plot(R2_Reis_Flow(:,1), R2_Reis_Flow(:,3), "x", 'LineWidth',2); hold off;
    title("Reis et al. Axyssimetric - Flow");
    legend("Levine and Schwinger","Reis et al.")
    xlabel("Frequency in Helmholtz number, ka");
    ylabel("|R|");
    xlim([0 1.8]);
    ylim([0.4 1.1]);
    set(gca,"FontSize",16);
    
%% 3.1 - Zhou Axisymmetric
nfig=nfig+1;
figure(nfig);
plot(LevineSchwinger(:,1), LevineSchwinger(:,3), 'LineWidth',2); hold on;
plot(R3_Zhou(:,1), R3_Zhou(:,3), "x", 'LineWidth',2); hold off;
    title("Zhou et al. Axyssimetric - No Flow");
    legend("Levine and Schwinger","Zhou et al.")
    xlabel("Frequency in Helmholtz number, ka");
    ylabel("|R|");
    xlim([0 1.8]);
    ylim([0.4 1.1]);
    set(gca,"FontSize",16);
    
%% 3.2 - Zhou Axisymmetric with FLOW
nfig=nfig+1;
figure(nfig);
plot(LevineSchwinger(:,1), LevineSchwinger(:,3), 'LineWidth',2); hold on;
plot(R3_Zhou_Flow(:,1), R3_Zhou_Flow(:,3), "x", 'LineWidth',2); hold off;
    title("Zhou et al. Axyssimetric - Flow");
    legend("Levine and Schwinger","Zhou et al.")
    xlabel("Frequency in Helmholtz number, ka");
    ylabel("|R|");
    xlim([0 1.8]);
    ylim([0.4 1.1]);
    set(gca,"FontSize",16);
    
%% All Results
nfig=nfig+1;
figure(nfig);
plot(LevineSchwinger(:,1), LevineSchwinger(:,3), "--", 'LineWidth',3); hold on;
plot(R1_2DChannel(:,1), R1_2DChannel(:,3), 'LineWidth',3);
plot(R1_2DChannel_Flow(:,1), R1_2DChannel_Flow(:,3), 'LineWidth',3);
plot(R2_Reis(:,1), R2_Reis(:,3), 'LineWidth',3);
plot(R2_Reis_Flow(:,1), R2_Reis_Flow(:,3), 'LineWidth',3);
plot(R3_Zhou(:,1), R3_Zhou(:,3), 'LineWidth',3);
plot(R3_Zhou_Flow(:,1), R3_Zhou_Flow(:,3), 'LineWidth',3); hold off;
    title("Reflection coefficient comparison");
    legend("Levine and Schwinger",...
        "2D Channel - No Flow",...
        "2D Channel - Flow",...
        "Reis et al. Axyssimetric - No Flow",...
        "Reis et al. Axyssimetric - Flow",...
        "Zhou - No Flow",...
        "Zhou - Flow");
    xlabel("Frequency in Helmholtz number, ka");
    ylabel("|R|");
    xlim([0 1.8]);
    ylim([0.3 1.1]);
    grid on;
    set(gca,"FontSize",16);
