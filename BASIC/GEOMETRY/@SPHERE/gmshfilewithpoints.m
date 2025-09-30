function varargout = gmshfilewithpoints(S,P,clS,clP,numbercenter,numberpoints,numberembeddedpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop,numbervolume,varargin)
% function G = gmshfilewithpoints(S,P,clS,clP,numbercenter,numberpoints,numberembeddedpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop,numbervolume)
% S : SPHERE
% P : POINT
% clS, clP : characteristic length

if ~iscell(P), P = {P}; end
if nargin<4 || isempty(clP), clP = clS; end
if isscalar(clP), clP = repmat(clP,1,length(P)); end

if nargin<5 || isempty(numbercenter), numbercenter = 1; end
if nargin<6 || isempty(numberpoints), numberpoints = numbercenter+(1:6); end
if nargin<7 || isempty(numberembeddedpoints), numberembeddedpoints = max([numbercenter,numberpoints])+(1:length(P)); end
if nargin<8 || isempty(numbercurves), numbercurves = 1:12; end
if nargin<9 || isempty(numbercurveloops), numbercurveloops = 1:8; end
if nargin<10 || isempty(numbersurfaces), numbersurfaces = 1:8; end
if nargin<11 || isempty(numbersurfaceloop), numbersurfaceloop = 1; end
if nargin<12 || isempty(numbervolume), numbervolume = 1; end

center = [S.cx,S.cy,S.cz];
PS = getvertices(S);

G = GMSHFILE();
G = createpoint(G,center,clS,numbercenter);
G = createpoints(G,PS,clS,numberpoints);
G = createspherecontour(G,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop);
G = createvolume(G,numbersurfaceloop,numbervolume);

G = createpoints(G,P,clP,numberembeddedpoints);
G = embedpointsinvolume(G,numberembeddedpoints,numbervolume);

if ischarin('recombine',varargin)
    G = recombinesurface(G);
end

varargout{1} = G;
varargout{2} = numbervolume;
