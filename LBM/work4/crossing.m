%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function was created to be used with the lattice Boltzmann code for 
% moving boundaries.
% Given the size of the main matrix, the direction of propagation and the
% two points that create a straight line, it finds the particles that
% crossed this line at any direction, so that:
%
% [d,c,l] = crossing(nl,nc,direc,xl,yl)
%
% Outputs:
% d is the vector containing the distances from the crossing points to the
% line segment defined by the ponts xl and yl.
% c and l are the column and line coordinates of the same points.
%
% Inputs:
% nl and nc are the number of lines and number of rows of the main matrix.
% direc is the direction of propagation acording to the D2Q9 model, defined 
% by the indexes 1->8:
%
%                       6   2   5
%                        \  |  /
%                       3 - . - 1
%                        /  |  \
%                       7   4   8
%
% xl and yl are the (1x2) vectors containg the coordinates of the initial
% and end points defining the boundary line segment.
%
% By Andrey R. da Silva                             McGill 06/2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s

function[desto,ce,le]=crossing(~,~,direc,xl,yl)



% direcao de propagacao 1 2 3 4 5 6 7 8
propx = [1 0 -1 0 1 -1 -1 1];
propy = [0 1 0 -1 1 1 -1 -1];
desto = [];
ce = [];
le = [];

% % definindo a matrix principal
% mat=ones(Nr,Mc);

for p =1:length(xl)-1
    %coeficientes da reta
    A = (yl(p+1)- yl(p))/(xl(p+1)-xl(p)+eps);
    B = yl(p+1)-A*xl(p+1);

    % definindo a linha paralela para o poligono
    xv=[xl(p:p+1) fliplr(xl(p:p+1)+propx(direc))];
    yv=[yl(p:p+1) fliplr(yl(p:p+1)+propy(direc))];


    % gera a caixa em torno do poligono
    x = floor(min(xv)):1:ceil(max(xv));
    y = floor(min(yv)):1:ceil(max(yv));

    % descreve os elementos dentro da caixa de forma vetorial
    tamanho=length(x)*length(y);
    c=zeros(1,tamanho);
    l=zeros(1,tamanho);
    fo=0;
    for u = 1 : length(x)
        for t = 1 : length(y)
             c(t+fo*length(y)) = x(1+fo);
             l(t+fo*length(y)) = y(t);
        end
        fo=fo+1;
    end

    % acha os pontos dentro do poligono
    [in, ~] = inpolygon(c,l,xv,yv);
    ji=find(in==1);
    c=c(ji);l=l(ji);
    norma=sqrt(2);
    % acha as distancias dos pontos que cruzaram a barreira da barreira em si,
    % para cada direcao

    if direc ==1 
        disto = c-(l-B)/(A+eps);
        elseif direc == 3
            disto = c-(l-B)/(A+eps);
        elseif direc == 2
            disto = l-(A*c+B);
        elseif direc == 4
            disto = l-(A*c+B);
        elseif direc == 5
            dx = c-(B-l+c)/(1-A);
            dy = l-(A*(l-c)-B)/(A-1);
            disto = sqrt(dx.^2+dy.^2)/norma;
        elseif direc == 7
            dx = c-(B-l+c)/(1-A);
            dy = l-(A*(l-c)-B)/(A-1);
            disto = sqrt(dx.^2+dy.^2)/norma;
        elseif direc == 8
            dx = c-(B-l-c)/(-A-1);
            dy = l-(A*l+A*c+B)/(1+A);
            disto = sqrt(dx.^2+dy.^2)/norma;
        elseif direc == 6
            dx = c-(B-l-c)/(-A-1);
            dy = l-(A*l+A*c+B)/(1+A);
            disto = sqrt(dx.^2+dy.^2)/norma;
        else
    end   
    desto=abs([desto disto]);
    ce = [ce c];
    le = [le l];
end
end
