function varargout = gmshfile(E,cl,numbercenter,numberpoints,numbercurves,numbercurveloop,numbersurface,varargin)
% function G = gmshfile(E,cl,numbercenter,numberpoints,numbercurves,numbercurveloop,numbersurface)
% E : ELLIPSE
% cl : characteristic length

if nargin<=2 || isempty(numbercenter), numbercenter = 1; end
if nargin<=3 || isempty(numberpoints), numberpoints = 2:5; end
if nargin<=4 || isempty(numbercurves), numbercurves = 1:4; end
if nargin<=5 || isempty(numbercurveloop), numbercurveloop = 1; numbersurface = 1; end
if nargin==6, numbersurface = []; end

semiaxes_lengths = [abs(E.a), abs(E.b)];
[~,maxidx] = max(semiaxes_lengths);

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
P = getvertices(E);
G = createpoint(G,[E.cx,E.cy,E.cz],cl,numbercenter);
G = createpoints(G,P,cl,numberpoints);
G = createellipsecontour(G,numbercenter,numberpoints,numbermajorpoints,numbercurves,numbercurveloop,varargin{:});
if ~isempty(numbersurface)
    G = createplanesurface(G,numbercurveloop,numbersurface);
    if ischarin('recombine',varargin)
        G = recombinesurface(G,numbersurface);
    end
end

varargout{1} = G;
varargout{2} = numbersurface;
