function [P,x,Xblock,mublock,I]=separatesamples2D(Xs,sepx,sepy,varargin)
% function [P,x,Xblock,mublock]=separatesamples(Xs,sepx,sepy)
% Xs: echantillons
% sepx: vecteur permettant la definition de n+1 intervalles selon la
% direction x
%     [-inf,z(1)] [z(1),z(2)] ... [z(end),Inf]
% sepy: idem pour la dimension y 
% P : vecteur contenant la probabilite de chaque intervalle
% x : vecteur contenant le decoupage associe de [0,1]
%       x = [0, P(1) , P(1)+P(2) , P(1)+P(2)+P(3) , ... , 1]
% Xblock : cellules contenant les echantillons de chacun des n+1 intervalles
% mublock : vecteur contenant la moyenne statistique de chaque block (esperance conditionnelle)

sepx = [-inf,sepx,inf];
sepy = [-inf,sepy,inf];
P = zeros(1,(length(sepx)-1)*(length(sepy)-1));
mublock = zeros(1,(length(sepx)-1)*(length(sepy)-1));
Xblock = cell(1,(length(sepx)-1)*(length(sepy)-1));

N = size(Xs,2);
I = zeros(size(Xs,2),1);


%rep=find(Xs(1,:)<=sepx(1) & Xs<=sepy(1));
%I(rep)=1;
%mublock(1)=mean(Xs(rep));
%Xblock{1} = Xs(rep);
%P(1) = numel(rep)/N;
for i=1:length(sepx)-1
    for j=1:length(sepy)-1
        num = (length(sepx)-1)*(j-1)+i;
        rep = find(Xs(1,:)>sepx(i) & Xs(1,:)<=sepx(i+1) & Xs(2,:)>sepy(j) & Xs(2,:)<=sepy(j+1));
        I(rep)=num;
        P(num) = numel(rep)/N;
        mublock(num)=mean(Xs(rep));
        Xblock{num} = Xs(rep);
    end
end

x = zeros(1,length(P)+1);
for i=1:length(P)
    x(i+1) = x(i)+P(i);
end    

if ischarin('display',varargin)
    figure(200)
    clf
    multipdfsampleplot(Xs,'npts',50)
    hold on
    for i=2:length(sepx)-1
        plot([sepx(i),sepx(i)],ylim,'r','LineWidth',3)
    end
    for i=2:length(sepy)-1
        plot(xlim,[sepy(i),sepy(i)],'r','LineWidth',3)
    end

end
