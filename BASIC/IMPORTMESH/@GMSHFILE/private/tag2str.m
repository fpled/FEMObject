function s = tag2str(tag)
% function s = tag2str(tag)
% Convert a tag to a char (character array) for Gmsh.
% tag: numeric/logical scalar (finite, real), scalar string (non-missing), or char row

if isnumeric(tag) || islogical(tag)
    if ~isscalar(tag)
        error('tag2str:NumericMustBeScalar', ...
              'Numeric tags must be scalar.');
    end
    if ~isreal(tag)
        error('tag2str:NumericMustBeReal', ...
              'Numeric tags must be real (no complex).');
    end
    if ~isfinite(tag)
        error('tag2str:NumericMustBeFinite', ...
              'Numeric tags must be finite (no NaN/Inf).');
    end
    s = num2str(tag);
    
elseif isstring(tag)
    if numel(tag)~=1
        error('tag2str:StringMustBeScalar', ...
              'String tags must be scalar.');
    end
    if ismissing(tag)
        error('tag2str:StringMissingNotAllowed', ...
              'Missing strings are not allowed as tags.');
    end
    s = char(tag);
    
elseif ischar(tag)
    if ~isrow(tag)
        error('tag2str:CharMustBeRow', ...
              'Char tags must be a 1-by-N row.');
    end
    s = tag;
    
else
    error('tag2str:UnsupportedElement', ...
          'Unsupported tag type "%s".', class(tag));
end

end
