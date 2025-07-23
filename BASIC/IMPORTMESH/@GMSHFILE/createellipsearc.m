function u = createellipsearc(u,numcenter,numpoints,numbermajorpoint,numbercurve)
% function u = createellipsearc(u,numcenter,numpoints,numbermajorpoint,numbercurve)

if length(numpoints)~=2
    error('An ellipse arc is defined by the 2 indices of the start and end points: numpoints = [start,end]')
end

if isempty(numbermajorpoint) || numpoints(1)==numbermajorpoint
    % Omit major axis point if the start point is a major axis point
    u = createentity(u,'Ellipse',[numpoints(1),numcenter,numpoints(2)],numbercurve);
else
    % Standard 4-point ellipse arc
    if ~isscalar(numbermajorpoint)
        error('An ellipse arc is defined by a scalar index of a major axis point: numbermajorpoint');
    end
    u = createentity(u,'Ellipse',[numpoints(1),numcenter,numbermajorpoint,numpoints(2)],numbercurve);
end
