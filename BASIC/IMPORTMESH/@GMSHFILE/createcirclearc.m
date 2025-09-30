function u = createcirclearc(u,numcenter,numpoints,numbercurve)
% function u = createcirclearc(u,numcenter,numpoints,numbercurve)

% if ~(isnumeric(numpoints) && numel(numpoints)==2)
%     error('A circle arc is defined by the 2 indices of the start and end points: numpoints = [start,end]')
% end
% 
% u = createentity(u,'Circle',[numpoints(1),numcenter,numpoints(2)],numbercurve);

% normalize the two point tags -> {startTag,endTag}
if isnumeric(numpoints) || islogical(numpoints)
    if numel(numpoints)~=2
        error('createcirclearc:PointsMustHaveTwo', ...
              'A circle arc requires exactly 2 point tags: [start, end].');
    end
    points = arrayfun(@tag2str, numpoints(:).', 'UniformOutput', false);

elseif isstring(numpoints)
    if numel(numpoints)~=2
        error('createcirclearc:PointsMustHaveTwo', ...
              'A circle arc requires exactly 2 point tags: [start, end].');
    end
    points = arrayfun(@tag2str, numpoints(:).', 'UniformOutput', false);

elseif ischar(numpoints)
    if size(numpoints,1)~=2
        error('createcirclearc:CharPointsMustBe2Rows', ...
              'Char array for points must have 2 rows (start and end).');
    end
    points = cellstr(numpoints).'; % {start,end}

elseif iscell(numpoints)
    if numel(numpoints)~=2
        error('createcirclearc:CellPointsMustHaveTwo', ...
              'Cell array for points must contain exactly 2 elements (start, end).');
    end
    points = cellfun(@tag2str, numpoints(:).', 'UniformOutput', false);

else
    error('createcirclearc:UnsupportedPointsType', ...
          'Unsupported type for "numpoints": %s', class(numpoints));
end

% center tag
center = tag2str(numcenter);

vals = [points(1), {center}, points(2)];

u = createentity(u,'Circle',vals,numbercurve);

end