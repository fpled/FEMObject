function s = CARTESIAN2D(varargin)
% function s = CARTESIAN2D(varargin)

if nargin==0
    s.name = 'Cartesian coodinates 2D';
    sp = SYSCOORD(2);
    s = class(s,'CARTESIAN2D',sp);
elseif nargin==1 && isa(varargin{1},'CARTESIAN2D')
    s = varargin{1};
else    
    s.name = 'Cartesian coodinates 2D';
    if isa(varargin{1},'VECTOR') && isa(varargin{2},'VECTOR')
        base = [getcompo(varargin{1}),getcompo(varargin{2})];
    end
    axis = {'X','Y'};
    sp = SYSCOORD(base,axis);
    s = class(s,'CARTESIAN2D',sp);
end
