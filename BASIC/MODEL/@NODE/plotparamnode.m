function varargout = plotparamnode(node,choix,varargin)
% function [Handles,leg] = plotparamnode(node,choix,varargin)
% choix : parametres à afficher
% Handles : handles vers les patchs d'elements
% leg : cellules contenant les legendes de chaque handle
%
% choix = 'nbddl'

paramplot = {};
Handles = [];
leg = {};

switch choix
    case 'nbddl'
        fun = @getnbddlpernode;
        
        
end

z = fun(node);
uz = unique(z);
nbgroup = length(uz);
for p=1:nbgroup
    groupnode{p}=find(z==uz(p));
    paramplot{p}=num2str(uz(p));
end


for p=1:nbgroup
    nodep = getnode(node,groupnode{p},'local');
    marker = getpointstyles(p);
    markersize=getcharin('MarkerSize',varargin,8);
    Htemp = plot(nodep,marker,varargin{:});
    Handles = [Handles,Htemp];
    leg = [leg,paramplot{p}];
end

axis off
axis image
if ~isempty(Handles)
    legend(Handles,leg{:});
end

if nargout>=1
    varargout{1}=Handles;
end
if nargout>=2
    varargout{2} = leg;
end


