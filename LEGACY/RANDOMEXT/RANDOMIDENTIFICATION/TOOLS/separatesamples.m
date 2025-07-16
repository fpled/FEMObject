function [P,x,Xblock,mublock,I]=separatesamples(Xs,sep,varargin)
% function [P,x,Xblock,mublock]=separatesamples(Xs,z)
% Xs: echantillons
% z: vecteur permettant la definition de n+1 intervalles
%     [-inf,z(1)] [z(1),z(2)] ... [z(end),Inf]
% P : vecteur contenant la probabilite de chaque intervalle
% x : vecteur contenant le decoupage associe de [0,1]
%       x = [0, P(1) , P(1)+P(2) , P(1)+P(2)+P(3) , ... , 1]
% Xblock : cellules contenant les echantillons de chacun des n+1 intervalles
% mublock : vecteur contenant la moyenne statistique de chaque block (esperance conditionnelle)

P = zeros(length(sep)+1,1);
mublock = zeros(length(sep)+1,1);
Xblock = cell(length(sep)+1,1);

N = numel(Xs);
I = zeros(numel(Xs),1);

rep=find(Xs<=sep(1));
I(rep)=1;
mublock(1)=mean(Xs(rep));
Xblock{1} = Xs(rep);
P(1) = numel(rep)/N;
for i=1:length(sep)-1
    rep = find(Xs>sep(i) & Xs<=sep(i+1));
    I(rep)=i+1;
%rep = find(Xs>=sep(i) & Xs<=sep(i+1));
    P(i+1) = numel(rep)/N;
    mublock(i+1)=mean(Xs(rep));
    Xblock{i+1} = Xs(rep);
end
rep = find(Xs>sep(end));
I(rep)=length(sep)+1;
%rep = find(Xs>=sep(end));
P(end) = numel(rep)/N;
mublock(end)=mean(Xs(rep));
Xblock{end} = Xs(rep);

x = zeros(1,length(P)+1);
for i=1:length(P)
    x(i+1) = x(i)+P(i);
end    

if ischarin('display',varargin)
    figure(200)
    clf
    pdfsampleplot(Xs,'b','npts',50,'bar')
    hold on
    for i=1:length(sep)
        plot([sep(i),sep(i)],ylim,'r','LineWidth',3)
    end
    if nargout==4
        for i=1:length(mublock)
            plot([mublock(i),mublock(i)],ylim,'g--','LineWidth',3)
        end
    end
end
