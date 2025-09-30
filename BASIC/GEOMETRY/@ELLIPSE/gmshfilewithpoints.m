function varargout = gmshfilewithpoints(E,P,clE,clP,numbercenter,numberpoints,numberembeddedpoints,numbercurves,numbercurveloop,numbersurface,varargin)
% function G = gmshfilewithpoints(E,P,clE,clP,numbercenter,numberpoints,numberembeddedpoints,numbercurves,numbercurveloop,numbersurface)
% E : ELLIPSE
% P : POINT
% clE, clP : characteristic length

if ~iscell(P), P = {P}; end
if nargin<4 || isempty(clP), clP = clE; end
if isscalar(clP), clP = repmat(clP,1,length(P)); end

if nargin<5 || isempty(numbercenter), numbercenter = 1; end
if nargin<6 || isempty(numberpoints), numberpoints = numbercenter+(1:4); end
if nargin<7 || isempty(numberembeddedpoints), numberembeddedpoints = max([numbercenter,numberpoints])+(1:length(P)); end
if nargin<8 || isempty(numbercurves), numbercurves = 1:4; end
if nargin<9 || isempty(numbercurveloop), numbercurveloop = 1; end
if nargin<10 || isempty(numbersurface), numbersurface = 1; end

center = [E.cx,E.cy,E.cz];
radii  = [E.a,E.b];
PE = getvertices(E);

[~,maxidx] = max(abs(radii));
switch maxidx
    case 1  % a is major (x-axis)
        numbermajorpoint = numberpoints(1); % -x
        % numbermajorpoint = numberpoints(3); % +x
    case 2  % b is major (y-axis)
        numbermajorpoint = numberpoints(2); % -y
        % numbermajorpoint = numberpoints(4); % +y
end
numbermajorpoints = repmat(numbermajorpoint,1,4); % Use the same major axis point for all 4 curves/arcs

G = GMSHFILE();
G = createpoint(G,center,clE,numbercenter);
G = createpoints(G,PE,clE,numberpoints);
G = createellipsecontour(G,numbercenter,numberpoints,numbermajorpoints,numbercurves,numbercurveloop);
G = createplanesurface(G,numbercurveloop,numbersurface);

G = createpoints(G,P,clP,numberembeddedpoints);
G = embedpointsinsurface(G,numberembeddedpoints,numbersurface);

if ischarin('recombine',varargin)
    G = recombinesurface(G,numbersurface);
end

varargout{1} = G;
varargout{2} = numbersurface;
