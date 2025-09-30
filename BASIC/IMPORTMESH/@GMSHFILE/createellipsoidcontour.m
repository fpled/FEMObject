function u = createellipsoidcontour(u,numbercenter,numberpoints,numbermajorpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop,varargin)
% function u = createellipsoidcontour(u,numcenter,numberpoints,numbermajorpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop)

if nargin<7, numbersurfaces = []; end
if nargin<8, numbersurfaceloop = []; end

% Naming convention for main points on the ellipsoid (order: +x, +y, +z, -x, -y, -z)
xp = numberpoints(1); % +x
yp = numberpoints(2); % +y
zp = numberpoints(3); % +z
xn = numberpoints(4); % -x
yn = numberpoints(5); % -y
zn = numberpoints(6); % -z

% 12 elliptical arcs (start point, end point, center, major axis point, arc number)
% The logic for major axis points should follow the non-colinear rule!
curves = [ ...
    zp, xp, numbercenter, numbermajorpoints(1),  numbercurves(1);  % +z, +x (north meridian)
    zp, yp, numbercenter, numbermajorpoints(2),  numbercurves(2);  % +z, +y
    zp, xn, numbercenter, numbermajorpoints(3),  numbercurves(3);  % +z, -x
    zp, yn, numbercenter, numbermajorpoints(4),  numbercurves(4);  % +z, -y
    
    xp, yp, numbercenter, numbermajorpoints(5),  numbercurves(5);  % +x, +y (equator)
    yp, xn, numbercenter, numbermajorpoints(6),  numbercurves(6);  % +y, -x
    xn, yn, numbercenter, numbermajorpoints(7),  numbercurves(7);  % -x, -y
    yn, xp, numbercenter, numbermajorpoints(8),  numbercurves(8);  % -y, +x
    
    zn, xp, numbercenter, numbermajorpoints(9),  numbercurves(9);  % -z, +x (south meridian)
    zn, yp, numbercenter, numbermajorpoints(10), numbercurves(10); % -z, +y
    zn, xn, numbercenter, numbermajorpoints(11), numbercurves(11); % -z, -x
    zn, yn, numbercenter, numbermajorpoints(12), numbercurves(12)  % -z, -y
];

for k=1:size(curves,1)
    u = createellipsearc(u,curves(k,3),[curves(k,1),curves(k,2)],curves(k,4),curves(k,5));
end

% 8 octants (3 arcs each)
octants = [ ...
    numbercurves(1),  numbercurves(5), -numbercurves(2);  % +z, +x, +y (north hemisphere)
    numbercurves(2),  numbercurves(6), -numbercurves(3);  % +z, +y, -x
    numbercurves(3),  numbercurves(7), -numbercurves(4);  % +z, -x, -y
    numbercurves(4),  numbercurves(8), -numbercurves(1);  % +z, -y, +x
    
    numbercurves(9),  numbercurves(5), -numbercurves(10); % -z, +x, +y (south hemisphere)
    numbercurves(10), numbercurves(6), -numbercurves(11); % -z, +y, -x
    numbercurves(11), numbercurves(7), -numbercurves(12); % -z, -x, -y
    numbercurves(12), numbercurves(8), -numbercurves(9)   % -z, -y, +x
];

for k=1:size(octants,1)
    u = createcurveloop(u,octants(k,:),numbercurveloops(k));
    if ~isempty(numbersurfaces)
        u = createsurface(u,numbercurveloops(k),numbersurfaces(k));
    end
end

if ~isempty(numbersurfaces) && ~isempty(numbersurfaceloop)
    u = createsurfaceloop(u,numbersurfaces,numbersurfaceloop);
end
