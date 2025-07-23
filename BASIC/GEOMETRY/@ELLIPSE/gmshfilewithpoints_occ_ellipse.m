function varargout = gmshfilewithpoints_occ_ellipse(E,P,clE,clP,numberembeddedpoints,numbercurve,numbercurveloop,numbersurface,varargin)
% function G = gmshfilewithpoints_occ_ellipse(E,P,clE,clP,numberembeddedpoints,numbercurve,numbercurveloop,numbersurface)
% E : ELLIPSE
% P : POINT
% clE, clP : characteristic length

if ~iscell(P), P = {P}; end
if nargin<=3 || isempty(clP), clP = clE; end
if isscalar(clP), clP = repmat(clP,1,length(P)); end

if nargin<=4 || isempty(numberembeddedpoints), numberembeddedpoints = 1+(1:length(P)); end
if nargin<=5 || isempty(numbercurve), numbercurve = 1; end
if nargin<=6 || isempty(numbercurveloop), numbercurveloop = 1; end
if nargin<=7 || isempty(numbersurface), numbersurface = 1; end

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
G = createellipse(G,[E.cx,E.cy,E.cz],[E.a,E.b],numbercurve);
G = setmeshsize(G,clE);
G = createcurveloop(G,numbercurve,numbercurveloop);
G = createplanesurface(G,numbercurveloop,numbersurface);

G = createpoints(G,P,clP,numberembeddedpoints);
G = embedpointsinsurface(G,numberembeddedpoints,numbersurface);

if ischarin('recombine',varargin)
    G = recombinesurface(G,numbersurface);
end

varargout{1} = G;
varargout{2} = numbersurface;
