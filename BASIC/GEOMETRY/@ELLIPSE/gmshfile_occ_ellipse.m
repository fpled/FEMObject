function varargout = gmshfile_occ_ellipse(E,cl,numbercurve,numbercurveloop,numbersurface,varargin)
% function G = gmshfile_occ_ellipse(E,cl,numbercurve,numbercurveloop,numbersurface)
% E : ELLIPSE
% cl : characteristic length

if nargin<=2 || isempty(numbercurve), numbercurve = 1; end
if nargin<=3 || isempty(numbercurveloop), numbercurveloop = 1; numbersurface = 1; end
if nargin==4, numbersurface = []; end

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
G = createellipse(G,[E.cx,E.cy,E.cz],[E.a,E.b],numbercurve);
G = setmeshsize(G,cl);
G = createcurveloop(G,numbercurve,numbercurveloop);
if ~isempty(numbersurface)
    G = createplanesurface(G,numbercurveloop,numbersurface);
    if ischarin('recombine',varargin)
        G = recombinesurface(G,numbersurface);
    end
end

varargout{1} = G;
varargout{2} = numbersurface;
