function u = createpoints(u,P,cl,numberpoints)
% function u = createpoints(u,P,cl,numberpoints)
%
if isa(P,'double')
    P = mat2cell(P,ones(1,size(P,1)),size(P,2));
end

npts = length(P);

if nargin < 3 || isempty(cl)
    cl = [];
elseif isscalar(cl)
    cl = repmat(cl,1,npts);
end

if nargin < 4 || isempty(numberpoints)
    numberpoints = [];
end

for k=1:npts
    Pk = P{k};
    if isempty(cl)
        if isempty(numberpoints)
            u = createpoint(u,Pk);
        else
            u = createpoint(u,Pk,[],numberpoints(k));
        end
    else
        if isempty(numberpoints)
            u = createpoint(u,Pk,cl(k));
        else
            u = createpoint(u,Pk,cl(k),numberpoints(k));
        end
    end
end
