function varargout = gmshfilewithpoints(D,P,clD,clP,numberpoints,numberembeddedpoints,numberlines,numberlineloop,numbersurface,varargin)
% function G = gmshfilewithpoints(D,P,clD,clP,numberpoints,numberembeddedpoints,numberlines,numberlineloop,numbersurface,numbersurfaceloop,numbervolume)
% D : DOMAIN
% P : POINT
% clD, clP : characteristic length

if ~iscell(P), P = {P}; end
if nargin<=3 || isempty(clP), clP = clD; end
if isscalar(clP), clP = repmat(clP,1,length(P)); end

if D.dim==2
    if nargin<=4 || isempty(numberpoints), numberpoints = 1:4; end
    if nargin<=5 || isempty(numberembeddedpoints), numberembeddedpoints = 4+(1:length(P)); end
    if nargin<=6 || isempty(numberlines), numberlines = 1:4; end
    if nargin<=7 || isempty(numberlineloop), numberlineloop = 1; end
    if nargin<=8 || isempty(numbersurface), numbersurface = 1; end
elseif D.dim==3
    if nargin<=4 || isempty(numberpoints), numberpoints = 1:8; end
    if nargin<=5 || isempty(numberembeddedpoints), numberembeddedpoints = 8+(1:length(P)); end
    if nargin<=6 || isempty(numberlines), numberlines = 1:12; end
    if nargin<=7 || isempty(numberlineloop), numberlineloop = 1:6; end
    if nargin<=8 || isempty(numbersurface), numbersurface = 1:6; end
    if nargin<=9 || isempty(varargin{1}), numbersurfaceloop = 1; end
    if nargin<=10 || isempty(varargin{2}), numbervolume = 1; end
end

G = GMSHFILE();
PD = getvertices(D);
if ischarin('extrude',varargin)
    G = createpoints(G,PD(1:4),clD,numberpoints(1:4));
else
    G = createpoints(G,PD,clD,numberpoints);
end
if D.dim==2
    G = createcontour(G,numberpoints,numberlines,numberlineloop);
    G = createplanesurface(G,numberlineloop,numbersurface);
elseif D.dim==3
    if ischarin('extrude',varargin)
        G = createcontour(G,numberpoints(1:4),numberlines(1:4),numberlineloop(1));
        vect = P{5}-P{1};
        [G,out] = extrude(G,vect,'Surface',numbersurface(1),varargin{:});
        numbertopsurface = [out,'[0]'];
        numbervolume     = [out,'[1]'];
        numberlateralsurfaces = {[out,'[2]'],[out,'[3]'],[out,'[4]'],[out,'[5]']};
        numbersurfaces = [{numbersurface}, {numbertopsurface}, numberlateralsurfaces{:}];
    else
        G = createdomaincontour(G,numberpoints,numberlines,numberlineloop,numbersurface,numbersurfaceloop);
        G = createvolume(G,numbersurfaceloop,numbervolume);
    end
end

G = createpoints(G,P,clP,numberembeddedpoints);
if D.dim==2
    G = embedpointsinsurface(G,numberembeddedpoints,numbersurface);
elseif D.dim==3
    G = embedpointsinvolume(G,numberembeddedpoints,numbervolume);  
end

if ischarin('recombine',varargin)
    if D.dim==2
        G = recombinesurface(G,numbersurface);
    elseif D.dim==3
        G = recombinesurface(G);
    end
end

varargout{1} = G;
if D.dim==2
    varargout{2} = numbersurface;
elseif D.dim==3
    varargout{2} = numbervolume;
    if ischarin('extrude',varargin)
        varargout{3} = numbersurfaces;
    end
end
