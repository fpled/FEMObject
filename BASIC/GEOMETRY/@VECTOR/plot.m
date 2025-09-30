function varargout = plot(u,varargin)
% function varargout = plot(u,varargin)

for l=1:nargin-1
    switch class(varargin{l})
        case 'POINT'
            v = varargin{l};
        case 'char'
            plotoptions = varargin{l};
        case 'double'
            ampl = varargin{l};
    end
end

if ~exist('plotoptions','var')
    plotoptions = 'b';
end
if ~exist('v','var')
    switch getindim(u)
        case 2
            v = POINT([0,0]);
        case 3
            v = POINT([0,0,0]);
        otherwise
            error('Wrong space dimension')
    end
end
if ~exist('ampl','var')
    ampl = 1;
end

u = double(u);
v = double(v);

hs = ishold;
hold on

switch length(v)
    case 2
        H = quiver(v(1),v(2),u(1),u(2),ampl,plotoptions);
    case 3
        H = quiver3(v(1),v(2),v(3),u(1),u(2),u(3),ampl,plotoptions);
    otherwise
        error('Wrong space dimension')
end

if ~hs
    hold off
end

if nargout
    varargout{1} = H;
end
