function u = recombinesurface(u,number,varargin)
% function u = recombinesurface(u,number[,angle])

if nargin<2 || isempty(number)
    number = ':';
end

options = '';
if nargin>2 && ~isempty(varargin)
    angle = varargin{1};
    options = [' = ' tag2str(angle)];
end

u = stepcounter(u);
s = ['Recombine Surface' tagsintobraces(number) options ';\n'];
u = addstring(u,s);
