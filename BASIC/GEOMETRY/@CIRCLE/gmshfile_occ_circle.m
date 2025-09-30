function varargout = gmshfile_occ_circle(C,cl,numbercurve,numbercurveloop,numbersurface,varargin)
% function G = gmshfile_occ_circle(C,cl,numbercurve,numbercurveloop,numbersurface)
% C : CIRCLE
% cl : characteristic length

if nargin<3 || isempty(numbercurve), numbercurve = 1; end
if nargin<4 || isempty(numbercurveloop), numbercurveloop = 1; end
if nargin<5, numbersurface = []; elseif isempty(numbersurface), numbersurface = 1; end

center = [C.cx,C.cy,C.cz];
radius = C.r;

v = [C.vx,C.vy];
n = [C.nx,C.ny,C.nz];
[dir,ang] = calcrotation_direction_angle(C,v,n);

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
G = createcircle(G,center,radius,numbercurve);
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
