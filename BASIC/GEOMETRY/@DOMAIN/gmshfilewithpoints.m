function varargout = gmshfilewithpoints(D,P,clD,clP,numberpoints,numberlines,numberlineloop,numberembeddedpoints,numbersurface,varargin)
% function G = gmshfilewithpoints(D,P,clD,clP,numberpoints,numberlines,numberlineloop,numberembeddedpoints,numbersurface)
% D : DOMAIN
% P : POINT
% clD, clP : characteristic length

if ~iscell(P)
    P = {P};
end
if nargin<=3
    clP = clD;
end
if length(clP)==1
    clP = repmat(clP,1,length(P));
end
if nargin<=4
    if getdim(D)==2
        numberpoints = 1:4;
    else
        numberpoints = 1:8;
    end
    numberlines = 1:4;
    numberlineloop = 5;
    if getdim(D)==2
        numberembeddedpoints = 4+(1:length(P));
    else
        numberembeddedpoints = 8+(1:length(P));
    end
    numbersurface = 1;
elseif nargin==8
    numbersurface = [];
end

G = GMSHFILE();
PD = getvertices(D);
G = createpoints(G,PD(1:4),clD,numberpoints);
G = createcontour(G,numberpoints(1:4),numberlines,numberlineloop);
if getdim(D)==3 || ~isempty(numbersurface)
    G = createplanesurface(G,numberlineloop,numbersurface);
    if getdim(D)==3
        vect = PD{5}-PD{1};
        G = extrude(G,vect,'Surface',numbersurface);
    end
end
if ~isempty(numbersurface)
    G = createpoints(G,P,clP,numberembeddedpoints);
    if getdim(D)==2
        G = embedpointsinsurface(G,numberembeddedpoints,numbersurface);
    elseif getdim(D)==3
        G = embedpointsinvolume(G,numberembeddedpoints,numbersurface);
    end
    if ischarin('recombine',varargin)
        G = recombinesurface(G,numbersurface);
    end
end

varargout{1} = G;
varargout{2} = numbersurface;
