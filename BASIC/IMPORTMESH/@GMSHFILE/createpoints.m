function u = createpoints(u,P,cl,numberpoints)
% function u = createpoints(u,P,cl,numberpoints)
%
if isa(P,'double')
    P = mat2cell(P,ones(1,size(P,1)),size(P,2));
end

if isscalar(cl)
    cl = repmat(cl,1,length(P));
end

for k=1:length(P)
    if isempty(cl)
        u = createpoint(u,P{k},[],numberpoints(k));
    else
        u = createpoint(u,P{k},cl(k),numberpoints(k));
    end
end
