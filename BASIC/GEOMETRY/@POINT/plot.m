function varargout = plot(u,varargin)
% function varargout = plot(u,varargin)

u = double(reshape2D(u));
u = [u,zeros(size(u,1),3-size(u,2))];
hold on

H = plot3(u(:,1),u(:,2),u(:,3),varargin{:});

if nargout
    varargout{1} = H;
end
