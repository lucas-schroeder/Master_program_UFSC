%% M?todo de Rayleigh-Ritz - Viga engastada-apoiada

%% Propriedades da viga
L=0.5;
b=0.02;
h=0.005;
Ar=b*h;
I=b*h^3/12;
pho=2700;
E=7.1e10;

%% Analise
n=[4 6 10]; %Numero de funcoes base em cada solucao

% Obs: Um numero de funcoes base maior que 10 apresenta problemas, o que
% pode ser devido a ma qualidade das funcoes base (elas passam a ter pouca
% diferenca entre elas para valores de i grandes).

x=0:0.01:L; %Pontos onde a resposta sera calculada

for m=1:length(n)
    N=n(m);
    clear M K d;
    
    for j=1:N
        for k=1:N
            % As constantes abaixo sao a norma das funcoes e sao usadas
            % para normatizar as funcoes (o que nem sempre eh necessario).
            C1=L^(j+5/2)*(1/(2*j+3)-2/(2*j+4)+1/(2*j+5))^0.5;
            C2=L^(k+5/2)*(1/(2*k+3)-2/(2*k+4)+1/(2*k+5))^0.5;
%             C1=1;
%             C2=1;
            % Matriz de massa e rigidez
            M(j,k)=pho*Ar*(L^(j+k+5)/(j+k+3)-2*L^(j+k+5)/(j+k+4)+L^(j+k+5)/(j+k+5))/(C1*C2);
            K(j,k)=E*I*(j*k*(j+1)*(k+1)*L^(j+k+1)/(j+k-1)-k*(j+2)*(j+1)*(k+1)*L^(j+k+1)/(j+k)-j*(j+1)*(k+2)*(k+1)*L^(j+k+1)/(j+k)+(j+2)*(j+1)*(k+2)*(k+1)*L^(j+k+1)/(j+k+1))/(C1*C2);
        end
    end

    [A,W]=eig(K,M);

    for j=1:N
        C1=L^(j+5/2)*(1/(2*j+3)-2/(2*j+4)+1/(2*j+5))^0.5;
        d(j,:)=x.^(j+1).*(L-x)/C1;
    end

    resultados(m).V=d'*A;
    resultados(m).fn=diag(W).^0.5/(2*pi);
    resultados(m).d=d;
end

% Funcoes base utilizadas na analise - 4 primeiras
d=resultados(1).d;
figure(1);
subplot(2,2,1);
plot(x,d(1,:));
axis([0 0.5 0 3]);
xlabel('x')
ylabel('d(x)')

subplot(2,2,2);
plot(x,d(2,:));
axis([0 0.5 0 3]);
xlabel('x')
ylabel('d(x)')

subplot(2,2,3);
plot(x,d(3,:));
axis([0 0.5 0 3]);
xlabel('x')
ylabel('d(x)')

subplot(2,2,4);
plot(x,d(4,:));
axis([0 0.5 0 3]);
xlabel('x')
ylabel('d(x)')



%Comparacao das frequencias naturais
for j=1:4
    for k=1:length(n)
        fn_N(j,k)=resultados(k).fn(j);
    end
end

figure(2);
semilogy(n,fn_N(1,:),'-o');
hold on;
semilogy(n,fn_N(2,:),'-dr');
semilogy(n,fn_N(3,:),'-sg');
semilogy(n,fn_N(4,:),'-<m');
xlabel('Numero de funcoes base')
ylabel('Freq. natural (Hz)')


%Comparacao das formas modais
figure(3);
subplot(2,2,1);
plot(x,resultados(1).V(:,1));
hold on;
plot(x,resultados(2).V(:,1),'--r');
plot(x,resultados(3).V(:,1),'-m');
axis([0 0.5 -5 5]);
xlabel('x [m]')
ylabel('Formal modal')
title('Primeiro Modo')

subplot(2,2,2);
plot(x,resultados(1).V(:,2));
hold on;
plot(x,resultados(2).V(:,2),'--r');
plot(x,resultados(3).V(:,2),'-m');
axis([0 0.5 -5 5]);
xlabel('x [m]')
ylabel('Formal modal')
title('Segundo Modo')

subplot(2,2,3);
plot(x,-resultados(1).V(:,3));
hold on;
plot(x,resultados(2).V(:,3),'--r');
plot(x,resultados(3).V(:,3),'-m');
axis([0 0.5 -5 5]);
xlabel('x [m]')
ylabel('Formal modal')
title('Terceiro Modo')

subplot(2,2,4);
plot(x,-resultados(1).V(:,4));
hold on;
plot(x,resultados(2).V(:,4),'--r');
plot(x,resultados(3).V(:,4),'-m');
axis([0 0.5 -7 7]);
xlabel('x [m]')
ylabel('Formal modal')
title('Quarto Modo')
legend('4 funcoes', '6 funcoes','10 funcoes')
