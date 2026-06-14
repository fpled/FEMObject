function u = createtransfinitevolume(u,number,varargin)
% function u = createtransfinitevolume(u,number[,values])

if nargin<2 || isempty(number)
    number = ':';
end

options = '';
if nargin>2 && ~isempty(varargin)
    values = varargin{1};
    options = [' = ' tagsintobraces(values)];
end

u = stepcounter(u);
s = ['Transfinite Volume ' tagsintobraces(number) options ';\n'];
u = addstring(u,s);
