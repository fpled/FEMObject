function varargout = gmshfile(S,cl,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop,numbervolume,varargin)
% function G = gmshfile(S,cl,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop,numbervolume)
% S : SPHERE
% cl : characteristic length

if nargin<=2 || isempty(numbercenter), numbercenter = 1; end
if nargin<=3 || isempty(numberpoints), numberpoints = 2:7; end
if nargin<=4 || isempty(numbercurves), numbercurves = 1:12; end
if nargin<=5 || isempty(numbercurveloops), numbercurveloops = 1:8; end
if nargin<=6 || isempty(numbersurfaces), numbersurfaces = 1:8; end
if nargin<=7 || isempty(numbersurfaceloop), numbersurfaceloop = 1; numbervolume = 1; end
if nargin==8, numbervolume = []; end

G = GMSHFILE();
P = getvertices(S);
G = createpoint(G,[S.cx,S.cy,S.cz],cl,numbercenter);
G = createpoints(G,P,cl,numberpoints);
G = createspherecontour(G,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop);
if ~isempty(numbervolume)
    G = createvolume(G,numbersurfaceloop,numbervolume);
end

if ischarin('recombine',varargin)
    G = recombinesurface(G);
end

varargout{1} = G;
varargout{2} = numbervolume;
