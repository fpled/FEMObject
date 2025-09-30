function u = createcylinderarccontour(u,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop,varargin)
% function u = createcylinderarccontour(u,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop)

if nargin<6, numbersurfaces = []; end
if nargin<7, numbersurfaceloop = []; end

n = numel(numberpoints)/2; % number of points per circle arc
if mod(numel(numberpoints),2)~=0
    error('createcylindercontour:PointCount', 'numberpoints length must be 2*n.');
end

% Points indices
numberbasepoints = numberpoints(1:n);     % n base points
numbertoppoints  = numberpoints(n+1:2*n); % n top points

% Curves (arcs/lines) indices
numberbasecurves = numbercurves(1:n+1);       % n-1 circle arcs + 2 radial lines at base
numbertopcurves  = numbercurves(n+2:2*n+2);   % n-1 circle arcs + 2 radial lines at top
numbervertlines  = numbercurves(2*n+3:3*n+2); % n side vertical lines between base and top 
numbercenterline = numbercurves(3*n+3);       % 1 center vertical line between base and top 

% Base circle arc (open, reverse orientation for outward normal)
u = createcirclearccontour(u,numbercenter(1),numberbasepoints,numberbasecurves,numbercurveloops(1),'reverse',-1);
if ~isempty(numbersurfaces)
    u = createplanesurface(u,numbercurveloops(1),numbersurfaces(1));
end

% Top circle arc (open)
u = createcirclearccontour(u,numbercenter(2),numbertoppoints,numbertopcurves,numbercurveloops(2));
if ~isempty(numbersurfaces)
    u = createplanesurface(u,numbercurveloops(2),numbersurfaces(2));
end

% Side vertical lines (n generatrices between base and top)
% seg_vert = [(1:n); (1:n)+n]';
% seg_vert = numberpoints(segvert);
seg_vert = [numberbasepoints(:), numbertoppoints(:)];
u = createlines(u,seg_vert,numbervertlines);

% Center vertical line (from base center to top center)
u = createline(u,[numbercenter(1),numbercenter(2)],numbercenterline);

% Lateral faces (n-1 quadrilaterals along the arc)
for i=1:(n-1)
    % Each lateral face: [base arc(i+1), vertical(i+1), -top arc(i+1), -vertical(i)]
    % curves = [ numbercurves(i+1), ...       % base arc (i+1 to i)
    %            numbercurves(2*n+2+i+1), ... % vertical up (i+1)
    %           -numbercurves(n+1+i+1), ...   % top arc (i to i+1)
    %           -numbercurves(2*n+2+i)];      % vertical down (i, reversed)
    curves = [ numberbasecurves(i+1), ... % base arc (i+1 to i) 
               numbervertlines(i+1), ...  % vertical up (i+1)
              -numbertopcurves(i+1), ...  % top arc (i to i+1, reversed)
              -numbervertlines(i)];       % vertical down (i, reversed)
    k = 2+i;
    u = createcurveloop(u,curves,numbercurveloops(k));
    if ~isempty(numbersurfaces)
        u = createsurface(u,numbercurveloops(k),numbersurfaces(k));
    end
end

% Radial faces (2 quadrilaterals)
% First (start) radial face (quadrilateral at start of arc)
% curves_first  = [ numbercurves(1), ...     % base radial (start)
%                   numbercurves(2*n+3), ... % side vertical (start)
%                  -numbercurves(n+2), ...   % top radial (start, reversed)
%                  -numbercurves(3*n+3) ];   % center vertical (reversed)
curves_first  = [ numberbasecurves(1), ... % base radial (start)
                  numbervertlines(1), ...  % side vertical (start)
                 -numbertopcurves(1), ...  % top radial (start, reversed)
                 -numbercenterline ];      % center vertical (reversed)
u = createcurveloop(u,curves_first,numbercurveloops(n+2));
if ~isempty(numbersurfaces)
    u = createplanesurface(u,numbercurveloops(n+2),numbersurfaces(n+2));
end

% Last (end) radial face (quadrilateral at end of arc)
% curves_last  = [ numbercurves(n+1), ...   % base radial (end)
%                  numbercurves(3*n+3), ... % center vertical
%                 -numbercurves(2*n+2), ... % top radial (end, reversed)
%                 -numbercurves(3*n+2) ];   % side vertical (end, reversed)
curves_last  = [ numberbasecurves(end), ... % base radial (end)
                 numbercenterline, ...      % center vertical
                -numbertopcurves(end), ...  % top radial (end, reversed)
                -numbervertlines(end) ];    % side vertical (end, reversed)
u = createcurveloop(u,curves_last,numbercurveloops(n+3));
if ~isempty(numbersurfaces)
    u = createplanesurface(u,numbercurveloops(n+3),numbersurfaces(n+3));
end

if ~isempty(numbersurfaces) && ~isempty(numbersurfaceloop)
    u = createsurfaceloop(u,numbersurfaces,numbersurfaceloop);
end
