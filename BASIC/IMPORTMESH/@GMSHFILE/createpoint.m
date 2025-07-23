function u = createpoint(u,P,cl,varargin)
% function u = createpoint(u,P,cl,number)

P = double(P);
if numel(P)==2
    P = [P,0];
end
if ~isempty(cl)
    P = [P,cl];
end
u = createentity(u,'Point',P,varargin{:});
