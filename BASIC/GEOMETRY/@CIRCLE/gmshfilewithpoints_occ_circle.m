function varargout = gmshfilewithpoints_occ_circle(C,P,clC,clP,numberembeddedpoints,numbercurve,numbercurveloop,numbersurface,varargin)
% function G = gmshfilewithpoints_occ_circle(C,P,clC,clP,numberembeddedpoints,numbercurve,numbercurveloop,numbersurface)
% C : CIRCLE
% P : POINT
% clC, clP : characteristic length

if ~iscell(P), P = {P}; end
if nargin<=3 || isempty(clP), clP = clC; end
if isscalar(clP), clP = repmat(clP,1,length(P)); end

if nargin<=4 || isempty(numberembeddedpoints), numberembeddedpoints = 1+(1:length(P)); end
if nargin<=5 || isempty(numbercurve), numbercurve = 1; end
if nargin<=6 || isempty(numbercurveloop), numbercurveloop = 1; end
if nargin<=7 || isempty(numbersurface), numbersurface = 1; end

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
G = createcircle(G,[C.cx,C.cy,C.cz],C.r,numbercurve);
G = setmeshsize(G,clC);
G = createcurveloop(G,numbercurve,numbercurveloop);
G = createplanesurface(G,numbercurveloop,numbersurface);

G = createpoints(G,P,clP,numberembeddedpoints);
G = embedpointsinsurface(G,numberembeddedpoints,numbersurface);

if ischarin('recombine',varargin)
    G = recombinesurface(G,numbersurface);
end

varargout{1} = G;
varargout{2} = numbersurface;
