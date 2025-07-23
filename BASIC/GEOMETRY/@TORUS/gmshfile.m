function varargout = gmshfile(T,cl,numbercenter,numberpoints,numbercurves,numbercurveloop,numbersurface,numbersurfaceloop,numbervolume,varargin)
% function G = gmshfile(T,cl,numbercenter,numberpoints,numbercurves,numbercurveloop,numbersurface,numbersurfaceloop,numbervolume)
% C : CYLINDER
% cl : characteristic length

if nargin<=2 || isempty(numbercenter), numbercenter = 1; end
if nargin<=3 || isempty(numberpoints), numberpoints = 2:7; end
if nargin<=4 || isempty(numbercurves), numbercurves = 1:12; end
if nargin<=5 || isempty(numbercurveloop), numbercurveloop = 1:8; end
if nargin<=6 || isempty(numbersurface), numbersurface = 1:8; end
if nargin<=7 || isempty(numbersurfaceloop), numbersurfaceloop = 1; numbervolume = 1; end
if nargin==8, numbervolume = []; end

G = GMSHFILE();
P = getvertices(T);
G = createpoint(G,[T.cx,T.cy,T.cz],cl,numbercenter);
G = createpoints(G,P,cl,numberpoints);
G = createtoruscontour(G,numbercenter,numberpoints,numbercurves,numbercurveloop,numbersurface);
if ~isempty(numbersurfaceloop)
    G = createsurfaceloop(G,numbersurface,numbersurfaceloop);
    if ~isempty(numbervolume)
        G = createvolume(G,numbersurfaceloop,numbervolume);
    end
end

if ischarin('recombine',varargin)
    G = recombinesurface(G);
end


varargout{1} = G;
varargout{2} = numbervolume;
