function varargout = gmshfile(E,cl,numbercenter,numberpoints,numbercurves,numbercurveloop,numbersurface,varargin)
% function G = gmshfile(E,cl,numbercenter,numberpoints,numbercurves,numbercurveloop,numbersurface)
% E : ELLIPSE
% cl : characteristic length

if nargin<3 || isempty(numbercenter), numbercenter = 1; end
if nargin<4 || isempty(numberpoints), numberpoints = numbercenter+(1:4); end
if nargin<5 || isempty(numbercurves), numbercurves = 1:4; end
if nargin<6 || isempty(numbercurveloop), numbercurveloop = 1; end
if nargin<7, numbersurface = []; elseif isempty(numbersurface), numbersurface = 1; end

center = [E.cx,E.cy,E.cz];
radii  = [E.a,E.b];
P = getvertices(E);

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
G = createpoint(G,center,cl,numbercenter);
G = createpoints(G,P,cl,numberpoints);
G = createellipsecontour(G,numbercenter,numberpoints,numbermajorpoints,numbercurves,numbercurveloop);
if ~isempty(numbersurface)
    G = createplanesurface(G,numbercurveloop,numbersurface);
    if ischarin('recombine',varargin)
        G = recombinesurface(G,numbersurface);
    end
end

varargout{1} = G;
varargout{2} = numbersurface;
