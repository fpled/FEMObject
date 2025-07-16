function varargout = plot(u,varargin)
% function varargout = plot(u,varargin)

numnode=ischarin('numnode',varargin);
varargin = delonlycharin('numnode',varargin);
if numnode
    color = getcharin('Color',varargin,'b');
    plotnumber(u,'Color',color);  
end
color = delcharin('Color',varargin);

H = plot(u.POINT,varargin{:});

if nargout>=1
    varargout{1}=H;
end


