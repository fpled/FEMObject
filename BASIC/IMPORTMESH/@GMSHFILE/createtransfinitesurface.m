function u = createtransfinitesurface(u,number,varargin)
% function u = createtransfinitesurface(u,number[,arrangement])
% function u = createtransfinitesurface(u,number[,values])
% function u = createtransfinitesurface(u,number[,values,arrangement])

if nargin<2 || isempty(number)
    number = ':';
end

options = '';
if nargin>2 && ~isempty(varargin)
    if isnumeric(varargin{1})
        values = varargin{1};
        options = [' = ' tagsintobraces(values)];
        if nargin>3 && ~isempty(varargin{2})
            arrangement = varargin{2};
            options = [' ' tag2str(arrangement)];
        end
    else
        arrangement = varargin{1};
        options = [' ' tag2str(arrangement)];
    end
end

u = stepcounter(u);
s = ['Transfinite Surface ' tagsintobraces(number) options ';\n'];
u = addstring(u,s);
