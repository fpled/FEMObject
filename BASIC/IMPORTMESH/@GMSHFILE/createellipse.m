function u = createellipse(u,center,radii,numbercurve,varargin)
% function u = createellipse(u,center,radii,numbercurve)
% function u = createellipse(u,center,radii,numbercurve,angles)

if length(center)~=3
    error('An ellipse is defined by the 3 coordinates of its center: center = [cx,cy,cz].')
end
if length(radii)~=2
    error('An ellipse is defined by its major (along the x-axis) and minor (along the y-axis) radii: radii = [r1,r2].')
end

if nargin==4 || isempty(varargin)
    u = createentity(u,'Ellipse',[center(:)' radii(:)'],numbercurve);
else
    angles = varargin{1};
    if length(angles)~=2
        error('An ellipse arc is defined by its start and end angles: angles = [angle1,angle2].')
    end
    u = createentity(u,'Ellipse',[center(:)' radii(:)' angles(:)'],numbercurve);
end
