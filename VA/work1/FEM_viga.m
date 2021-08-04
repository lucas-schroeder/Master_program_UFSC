%% M?todo de Elementos Finitos - Viga engastada-apoiada
%% Analise Modal
clear all;
close all;

%% Propriedades da viga
L=0.5;
b=0.02;
h=0.005;
A=b*h;
I=b*h^3/12;
pho=2700;
E=7.1e10;

%% Analise
n=[4 6 100]; %Numero de elementos da viga

for m=1:length(n)
    
    %Montando as matrizes do elemento
    a=L/2/n(m);
    Me=(pho*A*a/105)*[78 22*a 27 -13*a; 22*a 8*a^2 13*a -6*a^2; 27 13*a 78 -22*a; -13*a -6*a^2 -22*a 8*a^2];
    Ke=(E*I/(2*a^3))*[3 3*a -3 3*a; 3*a 4*a^2 -3*a 2*a^2; -3 -3*a 3 -3*a; 3*a 2*a^2 -3*a 4*a^2];
    
    %Montando as matrizes do sistem
    Nn=n(m)+1;
    Ngl=Nn*2;
    M=zeros(Ngl,Ngl);
    K=zeros(Ngl,Ngl);
    for j=1:n(m)
        Mee=zeros(Ngl,Ngl);
        Mee((2*j-1):(2*j+2),(2*j-1):(2*j+2))=Me;
        M=M+Mee;
        Kee=zeros(Ngl,Ngl);
        Kee((2*j-1):(2*j+2),(2*j-1):(2*j+2))=Ke;
        K=K+Kee;
    end
    
    %Aplicando as condicoes de contorno geometricas
    cc=[1 2 Ngl-1]; %Graus de Liberdade que devem ser restritos em ordem crescente.
    for j=1:length(cc)
        M(:,cc(j)-j+1) = []; % 1-1+1=1
        M(cc(j)-j+1,:) = [];
        K(:,cc(j)-j+1) = [];
        K(cc(j)-j+1,:) = [];
    end

    
    [Vc,W]=eig(K,M);
    
    resultados(m).fn=diag(W).^0.5/(2*pi);
    
    %Montando as formas modais
    %% Incluindo os GL das condicoes de contorno
    for j=1:length(cc)
        [s1,s2]=size(Vc);
        nline=zeros(1,s2);
        V1=Vc(1:(cc(j)-1),:);
        V2=Vc(cc(j):s1,:);
        Vc=cat(1,V1,nline,V2);
    end 
    %% Excluido as rotacoes
    V=[];
    for j=1:Ngl/2;
        V(j,:)=Vc(2*j-1,:);
    end
    
%     for j=1:Nn
%        V(1,j)=0;
%        V(2,j)=Vc(1,j);
%        V(3,j)=Vc(3,j);
%        V(4,j)=Vc(5,j);
%        V(5,j)=0;
%     end
    for j=1:size(W,1);
        if sum(V(:,j))>=0
            resultados(m).V(:,j)=V(:,j);
        else
            resultados(m).V(:,j)=-V(:,j);
        end
    end
end


%Comparacao das frequencias naturais
for j=1:4
    for k=1:length(n)
        fn_N(j,k)=resultados(k).fn(j);
    end
end

figure(1);
semilogy(n,fn_N(1,:),'-o');
hold on;
semilogy(n,fn_N(2,:),'-dr');
semilogy(n,fn_N(3,:),'-sg');
semilogy(n,fn_N(4,:),'-<m');
xlabel('Numero de funcoes base')
ylabel('Freq. natural (Hz)')


%Comparacao das formas modais
for j=1:length(n)
    resultados(j).x=0:L/n(j):L;
end
figure(3);
subplot(2,2,1);
plot(resultados(1).x,resultados(1).V(:,1));
hold on;
plot(resultados(2).x,resultados(2).V(:,1),'--r');
plot(resultados(3).x,resultados(3).V(:,1),'-k');
axis([0 0.5 -5 5]);
xlabel('x [m]')
ylabel('Formal modal')
title('Primeiro Modo')

subplot(2,2,2);
plot(resultados(1).x,resultados(1).V(:,2));
hold on;
plot(resultados(2).x,resultados(2).V(:,2),'--r');
plot(resultados(3).x,resultados(3).V(:,2),'-k');
axis([0 0.5 -5 5]);
xlabel('x [m]')
ylabel('Formal modal')
title('Segundo Modo')

subplot(2,2,3);
plot(resultados(1).x,resultados(1).V(:,3));
hold on;
plot(resultados(2).x,resultados(2).V(:,3),'--r');
plot(resultados(3).x,resultados(3).V(:,3),'-k');
axis([0 0.5 -5 5]);
xlabel('x [m]')
ylabel('Formal modal')
title('Terceiro Modo')

subplot(2,2,4);
plot(resultados(1).x,resultados(1).V(:,4));
hold on;
plot(resultados(2).x,resultados(2).V(:,4),'--r');
plot(resultados(3).x,resultados(3).V(:,4),'-k');
axis([0 0.5 -7 7]);
xlabel('x [m]')
ylabel('Formal modal')
title('Quarto Modo')
