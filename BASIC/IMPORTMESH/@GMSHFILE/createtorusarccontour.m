function u = createtorusarccontour(u,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop,varargin)
% function u = createtorusarccontour(u,numcenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop)

if nargin<6, numbersurfaces = []; end
if nargin<7, numbersurfaceloop = []; end

n = numel(numberpoints)/4; % number of minor-circle rings
if mod(numel(numberpoints),4)~=0
    error('createtorusarccontour:PointCount', 'numberpoints length must be 4*n.');
end

% Per-ring points (outer, top, inner, bottom) indices
numberpointO = @(i) numberpoints(4*(i-1)+1);
numberpointT = @(i) numberpoints(4*(i-1)+2);
numberpointI = @(i) numberpoints(4*(i-1)+3);
numberpointB = @(i) numberpoints(4*(i-1)+4);

% Curves (minor/major arcs) indices
numbercurve_minor = @(i,j) numbercurves(4*(i-1)+j); % j=1..4 (O->T, T->I, I->B, B->O)
numbercurve_majorO = @(i) numbercurves(4*n+i);      % i=1..n-1
numbercurve_majorT = @(i) numbercurves(4*n+1*(n-1)+i);
numbercurve_majorI = @(i) numbercurves(4*n+2*(n-1)+i);
numbercurve_majorB = @(i) numbercurves(4*n+3*(n-1)+i);

% Minor circle arcs on each ring: first 4*n curve tags
for i=1:n
    numbercenter_minor = numbercenter(3+i); % minor circle center
    numberpoints_minor = [numberpointO(i),numberpointT(i),numberpointI(i),numberpointB(i)];
    numbercurves_minor = numbercurves(4*(i-1)+(1:4)); % [O->T, T->I, I->B, B->O]
    % Create a curveloop for minor-circle ring 1 and ring n and cap them with plane surfaces
    if i==1
        u = createcirclecontour(u,numbercenter_minor,numberpoints_minor,numbercurves_minor,numbercurveloops(1),'reverse',-1);
        if ~isempty(numbersurfaces)
            u = createplanesurface(u,numbercurveloops(1),numbersurfaces(1));
        end
    elseif i==n
        u = createcirclecontour(u,numbercenter_minor,numberpoints_minor,numbercurves_minor,numbercurveloops(2));
        if ~isempty(numbersurfaces)
            u = createplanesurface(u,numbercurveloops(2),numbersurfaces(2));
        end
    else
        u = createcirclecontour(u,numbercenter_minor,numberpoints_minor,numbercurves_minor);
    end
end

% Major circle arcs: next 4*(n-1) curve tags, with no wrap
for i=1:n-1
    u = createcirclearc(u,numbercenter(1), [numberpointO(i),numberpointO(i+1)], numbercurve_majorO(i)); % outer circle arc
    u = createcirclearc(u,numbercenter(2), [numberpointT(i),numberpointT(i+1)], numbercurve_majorT(i)); % top circle arc
    u = createcirclearc(u,numbercenter(1), [numberpointI(i),numberpointI(i+1)], numbercurve_majorI(i)); % inner circle arc
    u = createcirclearc(u,numbercenter(3), [numberpointB(i),numberpointB(i+1)], numbercurve_majorB(i)); % bottom circle arc
end

% Lateral faces: 4 per ring (total 4*n)
for i=1:n-1
    % Quadrant 1: O->T
    curves = [ numbercurve_minor(i,1),   ... % (Oi -> Ti)
               numbercurve_majorT(i),    ... % (Ti -> Ti+1)
              -numbercurve_minor(i+1,1), ... % (Ti+1 -> Oi+1) reversed
              -numbercurve_majorO(i) ];      % (Oi+1 -> Oi) reversed
    k = 2+4*(i-1)+1;
    u = createcurveloop(u,curves,numbercurveloops(k));
    if ~isempty(numbersurfaces)
        u = createsurface(u,numbercurveloops(k),numbersurfaces(k));
    end

    % Quadrant 2: T->I
    curves = [ numbercurve_minor(i,2),   ...
               numbercurve_majorI(i),    ...
              -numbercurve_minor(i+1,2), ...
              -numbercurve_majorT(i) ];
    k = 2+4*(i-1)+2;
    u = createcurveloop(u,curves,numbercurveloops(k));
    if ~isempty(numbersurfaces)
        u = createsurface(u,numbercurveloops(k),numbersurfaces(k));
    end

    % Quadrant 3: I->B
    curves = [ numbercurve_minor(i,3),   ...
               numbercurve_majorB(i),    ...
              -numbercurve_minor(i+1,3), ...
              -numbercurve_majorI(i) ];
    k = 2+4*(i-1)+3;
    u = createcurveloop(u,curves,numbercurveloops(k));
    if ~isempty(numbersurfaces)
        u = createsurface(u,numbercurveloops(k),numbersurfaces(k));
    end
    
    % Quadrant 4: B->O
    curves = [ numbercurve_minor(i,4),   ...
               numbercurve_majorO(i),    ...
              -numbercurve_minor(i+1,4), ...
              -numbercurve_majorB(i) ];
    k = 2+4*(i-1)+4;
    u = createcurveloop(u,curves,numbercurveloops(k));
    if ~isempty(numbersurfaces)
        u = createsurface(u,numbercurveloops(k),numbersurfaces(k));
    end
end

if ~isempty(numbersurfaces) && ~isempty(numbersurfaceloop)
    u = createsurfaceloop(u, numbersurfaces, numbersurfaceloop);
end

end
