function u = createellipsearc(u,numcenter,numpoints,numbermajorpoint,numbercurve)
% function u = createellipsearc(u,numcenter,numpoints,numbermajorpoint,numbercurve)

% if ~(isnumeric(numpoints) && numel(numpoints)==2)
%     error('An ellipse arc is defined by the 2 indices of the start and end points: numpoints = [start,end]')
% end
% 
% if isempty(numbermajorpoint) || numpoints(1)==numbermajorpoint
%     % Omit major axis point if the start point is a major axis point
%     u = createentity(u,'Ellipse',[numpoints(1),numcenter,numpoints(2)],numbercurve);
% else
%     % Standard 4-point ellipse arc
%     if ~isscalar(numbermajorpoint)
%         error('An ellipse arc is defined by a scalar index of a major axis point: numbermajorpoint');
%     end
%     u = createentity(u,'Ellipse',[numpoints(1),numcenter,numbermajorpoint,numpoints(2)],numbercurve);
% end


% normalize the two point tags -> {startTag,endTag}
if isnumeric(numpoints) || islogical(numpoints)
    if numel(numpoints)~=2
        error('createellipsearc:PointsMustHaveTwo', ...
              'An ellipse arc requires exactly 2 point tags: [start, end].');
    end
    points = arrayfun(@tag2str, numpoints(:).', 'UniformOutput', false);

elseif isstring(numpoints)
    if numel(numpoints)~=2
        error('createellipsearc:PointsMustHaveTwo', ...
              'An ellipse arc requires exactly 2 point tags: [start, end].');
    end
    points = arrayfun(@tag2str, numpoints(:).', 'UniformOutput', false);

elseif ischar(numpoints)
    if size(numpoints,1)~=2
        error('createellipsearc:CharPointsMustBe2Rows', ...
              'Char array for points must have 2 rows (start and end).');
    end
    points = cellstr(numpoints).'; % {start,end}

elseif iscell(numpoints)
    if numel(numpoints)~=2
        error('createellipsearc:CellPointsMustHaveTwo', ...
              'Cell array for points must contain exactly 2 elements (start, end).');
    end
    points = cellfun(@tag2str, numpoints(:).', 'UniformOutput', false);

else
    error('createellipsearc:UnsupportedPointsType', ...
          'Unsupported type for "numpoints": %s', class(numpoints));
end

% center tag
center = tag2str(numcenter);

% If major axis point equals start point, Gmsh allows 3-arg form
if nargin<4 || isempty(numbermajorpoint)
    vals = [points(1), {center}, points(2)]; % 3-arg form
else
    % major axis point tag
    majorpoint = tag2str(numbermajorpoint);
    if strcmp(majorpoint, points{1})
        vals = [points(1), {center}, points(2)]; % omit major axis point if equal to start point
    else
        vals = [points(1), {center, majorpoint}, points(2)]; % full 4-arg form
    end
end

u = createentity(u,'Ellipse',vals,numbercurve);

end
