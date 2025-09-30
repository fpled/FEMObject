function varargout = gmshfile(S,cl,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop,numbervolume,varargin)
% function G = gmshfile(S,cl,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop,numbervolume)
% S : SPHERE
% cl : characteristic length

if nargin<3 || isempty(numbercenter), numbercenter = 1; end
if nargin<4 || isempty(numberpoints), numberpoints = numbercenter+(1:6); end
if nargin<5 || isempty(numbercurves), numbercurves = 1:12; end
if nargin<6 || isempty(numbercurveloops), numbercurveloops = 1:8; end
if nargin<7 || isempty(numbersurfaces), numbersurfaces = 1:8; end
if nargin<8 || isempty(numbersurfaceloop), numbersurfaceloop = 1; end
if nargin<9, numbervolume = []; elseif isempty(numbervolume), numbervolume = 1; end

center = [S.cx,S.cy,S.cz];
P = getvertices(S);

G = GMSHFILE();
G = createpoint(G,center,cl,numbercenter);
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
