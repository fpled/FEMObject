function u = createcylindercontour(u,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop)
% function u = createcylindercontour(u,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop)

n = length(numberpoints)/2;  % number of points per circle

% Points indices
numberbasepoints = numberpoints(1:n);     % n base points
numbertoppoints  = numberpoints(n+1:2*n); % n top points

% Curves (arcs/lines) indices
numberbasecurves = numbercurves(1:n);       % n base circle arcs
numbertopcurves  = numbercurves(n+1:2*n);   % n top circle arcs
numbervertlines  = numbercurves(2*n+1:3*n); % n vertical lines between base and top

% Base circle (closed, reverse orientation for outward normal)
u = createcirclecontour(u,numbercenter(1),numberbasepoints,numberbasecurves,numbercurveloops(1),'reverse',-1);
u = createplanesurface(u,numbercurveloops(1),numbersurfaces(1));

% Top circle (closed)
u = createcirclecontour(u,numbercenter(2),numbertoppoints,numbertopcurves,numbercurveloops(2));
u = createplanesurface(u,numbercurveloops(2),numbersurfaces(2));

% Vertical lines (n generatrices between base and top)
% segvert = [(1:n); (1:n)+n]';
% segvert = numberpoints(segvert);
segvert = [numberbasepoints(:), numbertoppoints(:)];
u = createlines(u,segvert,numbervertlines);

% Lateral faces (n quadrilaterals along the circle)
for i=1:n
    j = mod(i,n) + 1; % wraps to 1 when i = n, ensuring the circle/loop closes correctly
    % Each lateral face: [base arc(i), vertical(j), -top arc(i), -vertical(i)]
    % curves = [ numbercurves(i), ...     % base arc (j to i)
    %            numbercurves(2*n+j), ... % vertical up (j)
    %           -numbercurves(n+i), ...   % top arc (i to j, reversed)
    %           -numbercurves(2*n+i)];    % vertical down (i, reversed)
    curves = [ numberbasecurves(i), ... % base arc (j to i)
               numbervertlines(j), ...  % vertical up (j)
              -numbertopcurves(i), ...  % top arc (i to j, reversed)
              -numbervertlines(i)];     % vertical down (i, reversed)
    k = 2+i;
    u = createcurveloop(u,curves,numbercurveloops(k));
    u = createsurface(u,numbercurveloops(k),numbersurfaces(k));
end
u = createsurfaceloop(u,numbersurfaces,numbersurfaceloop);
