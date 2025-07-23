function s = SYSCOORD(varargin)
% function s = SYSCOORD(varargin)

if nargin==0
    s = SYSCOORD(2);
elseif nargin==1 && isa(varargin{1},'char')
    switch varargin{1}
        case 'UNID'
            dim = 1;
        case 'PLANE'
            dim = 2;
        case 'TRID'
            dim = 3;
    end
    axis = {'X','Y','Z'};
    base = MYDOUBLEND(eye(dim));
    
    s.dim = dim;
    s.indim = dim;
    s.axis = axis(1:dim);
    
    s = class(s,'SYSCOORD',base);
    
elseif nargin==1 && isa(varargin{1},'double') && numel(varargin{1})==1
    dim = varargin{1};
    axis = {'X','Y','Z'};
    axis = axis(1:dim);
    base = MYDOUBLEND(eye(dim));
    
    s.dim = dim;
    s.indim = dim;
    s.axis = axis;
    
    s = class(s,'SYSCOORD',base);
    
elseif nargin==2
    
    base = MYDOUBLEND(varargin{1});
    axis = varargin{2};
    
    s.dim = size(base,2);
    s.indim = size(base,1);
    s.axis = axis;
    
    s = class(s,'SYSCOORD',base);
    
else
    error('Wrong input arguments')
end
