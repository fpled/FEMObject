function varargout = gmshfilewithpoints(E,P,clE,clP,numbercenter,numberpoints,numberembeddedpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop,numbervolume,varargin)
% function G = gmshfilewithpoints(E,P,clE,clP,numbercenter,numberpoints,numberembeddedpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop,numbervolume)
% E : ELLIPSOID
% P : POINT
% clE, clP : characteristic length

if ~iscell(P), P = {P}; end
if nargin<=3 || isempty(clP), clP = clE; end
if isscalar(clP), clP = repmat(clP,1,length(P)); end

if nargin<=4 || isempty(numbercenter), numbercenter = 1; end
if nargin<=5 || isempty(numberpoints), numberpoints = 2:7; end
if nargin<=6 || isempty(numberembeddedpoints), numberembeddedpoints = numberpoints(end)+(1:length(P)); end
if nargin<=7 || isempty(numbercurves), numbercurves = 1:12; end
if nargin<=8 || isempty(numbercurveloops), numbercurveloops = 1:8; end
if nargin<=9 || isempty(numbersurfaces), numbersurfaces = 1:8; end
if nargin<=10 || isempty(numbersurfaceloop), numbersurfaceloop = 1; end
if nargin<=11 || isempty(numbervolume), numbervolume = 1; end

semiaxes_lengths = [abs(E.a), abs(E.b), abs(E.c)];
[~,maxidx] = max(semiaxes_lengths);

switch maxidx
    case 1  % a is major (x-axis)
        numbermajorpoint = numberpoints(1); % -x
        % numbermajorpoint = numberpoints(4); % +x
    case 2  % b is major (y-axis)
        numbermajorpoint = numberpoints(2); % -y
        % numbermajorpoint = numberpoints(5); % +y
    case 3  % c is major (z-axis)
        numbermajorpoint = numberpoints(3); % -z
        % numbermajorpoint = numberpoints(6); % +z
end
numbermajorpoints = repmat(numbermajorpoint,1,12); % Use the same major axis point for all 12  curves/arcs

G = GMSHFILE();
PE = getvertices(E);
G = createpoint(G,[E.cx,E.cy,E.cz],clE,numbercenter);
G = createpoints(G,PE,clE,numberpoints);
G = createellipsoidcontour(G,numbercenter,numberpoints,numbermajorpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop);
G = createvolume(G,numbersurfaceloop,numbervolume);

G = createpoints(G,P,clP,numberembeddedpoints);
G = embedpointsinvolume(G,numberembeddedpoints,numbervolume);

if ischarin('recombine',varargin)
    G = recombinesurface(G);
end

varargout{1} = G;
varargout{2} = numbervolume;
