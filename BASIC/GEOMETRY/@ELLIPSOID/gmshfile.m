function varargout = gmshfile(E,cl,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop,numbervolume,varargin)
% function G = gmshfile(E,cl,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop,numbervolume)
% E : ELLIPSOID
% cl : characteristic length

if nargin<3 || isempty(numbercenter), numbercenter = 1; end
if nargin<4 || isempty(numberpoints), numberpoints = numbercenter+(1:6); end
if nargin<5 || isempty(numbercurves), numbercurves = 1:12; end
if nargin<6 || isempty(numbercurveloops), numbercurveloops = 1:8; end
if nargin<7 || isempty(numbersurfaces), numbersurfaces = 1:8; end
if nargin<8 || isempty(numbersurfaceloop), numbersurfaceloop = 1; end
if nargin<9, numbervolume = []; elseif isempty(numbervolume), numbervolume = 1; end

center = [E.cx,E.cy,E.cz];
P = getvertices(E);

semiaxes_lengths = [abs(E.a),abs(E.b),abs(E.c)];
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
G = createpoint(G,center,cl,numbercenter);
G = createpoints(G,P,cl,numberpoints);
G = createellipsoidcontour(G,numbercenter,numberpoints,numbermajorpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop);
if ~isempty(numbervolume)
    G = createvolume(G,numbersurfaceloop,numbervolume);
end

if ischarin('recombine',varargin)
    G = recombinesurface(G);
end

varargout{1} = G;
varargout{2} = numbervolume;
