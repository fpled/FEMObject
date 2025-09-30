function u = setfieldattribute(u,name,value,number)
% function u = setfieldattribute(u,name,value,number)

if nargin<4 || isempty(number)
    number = 1;
end

u = stepcounter(u);

if ischar(value)
    val = value;
elseif isstring(value) && isscalar(value)
    val = ['"' char(value) '"'];
elseif (isnumeric(value) || islogical(value)) && isscalar(value)
    val = num2str(value);
elseif (isnumeric(value) || islogical(value)) && ~isscalar(value)
    val = tagsintobraces(value);
else
    error('setfieldattribute:UnsupportedType', ...
        'Unsupported value type for setfieldattribute.');
end

s = ['Field[' num2str(number) '].' name ' = ' val ';\n'];
u = addstring(u,s);
