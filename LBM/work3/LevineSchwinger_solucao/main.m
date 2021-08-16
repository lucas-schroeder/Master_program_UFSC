%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% finds the reflectance using the Levine and Schwinger solution %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc



dens =1;
c =340;
raio = 0.025/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do explicit magnitude calculation using Eq. (VI.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = linspace(0,3.8,300); % selects the number of points on the curve

sum1 = zeros(size(z));
for i=1:length(z),
 	sum1(i) = quadl('int3',eps,20,[],[],z(i));
end
RR = (pi*z).^(0.5).*exp(-z+(sum1./pi));
RR(1) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do explicit end correction calculation using Eq. (VI.4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum1 = size(z);
sum2 = size(z);
sum1(1) = 0;
for i=2:length(z),
  sum1(i) = quadl('int4',eps,z(i)-0.0000001,[],[],z(i));
end

sum2(1) = 0;
for i=2:length(z),
  sum2(i) = quadl('int2',eps,600,[],[],z(i));
end
loa = (sum1+sum2)./pi;
loa(1) = 0.6133;          %% loa = l/a where a is the radius of the pipe

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the complex R according to Eq. (VI.2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rc = -RR.*exp(2*j*z.*loa);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tranformin the complex reflectance into complex impedance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Z = dens*c*(1+Rc)./(1-Rc);
Z = Z/(dens*c);



%save_ZL_norm Z z

figure(1)
plot(z,real(Z/(dens*c)),'r')
hold on
plot(z,abs(imag(Z)/(dens*c)),'b')
title('Radiation impedance - Levine and Schwinger exact solution for unflanged pipe ')
legend('Z real part','Z imaginary part')
xlabel('Helmholtz number ka')
ylabel('Impedance Z')


figure(2)
plot(z,RR)
title('Reflection coefficient at the open end - Levine and Schwinger exact solution for an unflanged pipe ')
xlabel('Helmholtz number ka')
ylabel('Magnitude of R')
%% Save results
LevineSchwinger = [z.' Z.' RR.'];
save('LevineSchwinger.mat','LevineSchwinger');
