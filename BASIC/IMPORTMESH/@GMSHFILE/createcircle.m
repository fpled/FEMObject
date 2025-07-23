function u = createcircle(u,center,radius,numbercurve,varargin)
% function u = createcircle(u,center,radius,numbercurve)
% function u = createcircle(u,center,radius,numbercurve,angles)

if length(center)~=3
    error('A circle is defined by the 3 coordinates of its center: center = [cx,cy,cz].')
end
if isempty(radius) || ~isscalar(radius)
    error('A circle is defined by its scalar radius.');
end

if nargin==4 || isempty(varargin)
    u = createentity(u,'Circle',[center(:)' radius],numbercurve);
else
    angles = varargin{1};
    if length(angles)~=2
        error('A circle arc is defined by its start and end angles: angles = [angle1,angle2].')
    end
    u = createentity(u,'Circle',[center(:)' radius angles(:)'],numbercurve);
end
