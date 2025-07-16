function plotVwithname(u,scanm,name,varargin)
% function plotVwithname(u,scanm,name,varargin)
% affichage des modes de la PCRADIALMATRIX
% scanm : ensembles des modes à tracer
% on utilise plot(V{i},varargin{:})
% varargin : arguments de la fonction plot
if isempty(scanm) | nargin==1
    scanm = 1:u.m;
end
m = length(scanm);
nl=floor(sqrt(m));nc=ceil(m/nl);
for i=1:m
fullsubplot(nl,nc,i);

V = u.V{scanm(i)};

plot(V,varargin{:});
axis off
title([ name '_{' num2str(scanm(i)) '}'],'FontSize',16)
end
