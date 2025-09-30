function u = createdomaincontour(u,numberpoints,numberlines,numbercurveloops,varargin)
% function u = createdomaincontour(u,numberpoints,numberlines,numbercurveloops)
% function u = createdomaincontour(u,numberpoints,numberlines,numbercurveloops,numbersurfaces,numbersurfaceloop)

reverse = getcharin('reverse',varargin,1);
assert(ismember(reverse,[-1,1]),'''reverse'' must be +1 or -1');

n = length(numberpoints);

if n==4
    % 2D quadrilateral (rectangle): 4 points, 4 lines, 1 line loop
    % seg = [1:n; 2:n,1]';
    % seg = [1 2 3 4; 2 3 4 1]';
    % seg = [1 2; 2 3; 3 4; 4 1];
    % seg = numberpoints(seg);
    % u = createlines(u,seg,numberlines);
    % if ~isempty(numbercurveloops)
    %     numberlines = reverse * numberlines;
    %     u = createcurveloop(u,numberlines,numbercurveloops);
    % end
    u = createcontour(u,numberpoints,numberlines,numbercurveloops,'reverse',reverse);
elseif n==8
    if nargin<5, numbersurfaces=[]; else, numbersurfaces = varargin{1}; end
    if nargin<6, numbersurfaceloop=[]; else, numbersurfaceloop = varargin{2}; end

    % 3D hexahedron (box): 8 points, 12 lines, 6 line loops, 6 plane surfaces, 1 surface loop
    % 8 vertices (points): 1-4 bottom, 5-8 top
    % 12 lines (edges)   : 1-4 bottom, 5-8 top, 9-12 verticals
    % bottom = 1:4; top = 5:8;
    % seg = [bottom, top, bottom; ...
    %        circshift(bottom,-1), circshift(top,-1), top]';
    % seg = [1 2 3 4  5 6 7 8  1 2 3 4;
    %        2 3 4 1  6 7 8 5  5 6 7 8]';
    seg = [1 2; 2 3; 3 4; 4 1; ... % bottom
           5 6; 6 7; 7 8; 8 5; ... % top
           1 5; 2 6; 3 7; 4 8];    % verticals
    seg = numberpoints(seg);
    u = createlines(u,seg,numberlines);
    
    % 6 faces (signed for orientation)
    faces = [1 2 3 4;     % bottom (1-2-3-4)
             5 6 7 8;     % top    (5-6-7-8)
             1 10 -5 -9;  % front  (1-2-6-5)
             2 11 -6 -10; % left   (2-3-7-6)
             3 12 -7 -11; % back   (3-4-8-7)
             4 9 -8 -12]; % right  (4-1-5-8)
    for k=1:size(faces,1)
        lines = numberlines(abs(faces(k,:)));
        orient = sign(faces(k,:));
        lines = lines .* orient;
        u = createlineloop(u,lines,numbercurveloops(k));
        if ~isempty(numbersurfaces)
            u = createplanesurface(u,numbercurveloops(k),numbersurfaces(k));
        end
    end
    
    if ~isempty(numbersurfaces) && ~isempty(numbersurfaceloop)
        numbersurfaces = reverse * numbersurfaces;
        u = createsurfaceloop(u,numbersurfaces,numbersurfaceloop);
    end
else
    error('createdomaincontour: Only supports 4 or 8 points.');
end
