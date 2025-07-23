function u = createdisk(u,center,radius,numbersurface,varargin)
% function u = createdisk(u,center,radius,numbersurface)
% function u = createdisk(u,center,radius,numbersurface,radiusy)

if length(center)~=3
    error('A disk is defined by the 3 coordinates of its center: center = [cx,cy,cz].')
end
if isempty(radius) || ~isscalar(radius)
    error('A disk is defined by its scalar radius.');
end

u = createentity(u,'Disk',[center(:)' radius varargin{:}],numbersurface);
