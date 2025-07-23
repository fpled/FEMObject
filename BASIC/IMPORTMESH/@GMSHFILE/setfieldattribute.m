function u = setfieldattribute(u,name,value,number)
% function u = setfieldattribute(u,name,value,number)

if nargin<4 || isempty(number)
    number = 1;
end

u = stepcounter(u);

if ischar(value) || (isstring(value) && isscalar(value))
    valstr = ['"' char(value) '"'];
elseif isnumeric(value) && numel(value)>1
    valstr = valuesintobraces(value);
elseif isnumeric(value) && isscalar(value)
    valstr = num2str(value);
else
    error('Unsupported value type for setfieldattribute.');
end

s = ['Field[' num2str(number) '].' name ' = ' valstr ' ;\n'];
u = addstring(u,s);
