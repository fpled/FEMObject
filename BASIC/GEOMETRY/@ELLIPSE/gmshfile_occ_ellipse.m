function varargout = gmshfile_occ_ellipse(E,cl,numbercurve,numbercurveloop,numbersurface,varargin)
% function G = gmshfile_occ_ellipse(E,cl,numbercurve,numbercurveloop,numbersurface)
% E : ELLIPSE
% cl : characteristic length

if nargin<3 || isempty(numbercurve), numbercurve = 1; end
if nargin<4 || isempty(numbercurveloop), numbercurveloop = 1; end
if nargin<5, numbersurface = []; elseif isempty(numbersurface), numbersurface = 1; end

center = [E.cx,E.cy,E.cz];
radii  = [E.a,E.b];

v = [E.vx,E.vy];
n = [E.nx,E.ny,E.nz];
[dir,ang] = calcrotation_direction_angle(E,v,n);

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
G = createellipse(G,center,radii,numbercurve);
if abs(ang) > eps
    G = rotate(G,dir,center,ang,'Curve',numbercurve);
end
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
