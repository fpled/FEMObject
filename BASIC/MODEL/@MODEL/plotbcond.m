function varargout = plotbcond(M,varargin)
% function varargout = plotbcond(M,varargin)

BC = getbc(M.BCOND);
paramplot = cell(1,length(BC));
groupnode = cell(1,length(BC));

for i=1:length(BC)
    if BC{i}.type==0
        name = 'Dirichlet';
        [numnode,ddl] = findnodeinnumddl(M,BC{i}.ddlbloque);
        groupnode{i} = numnode;
    elseif BC{i}.type==1
        name = 'Periodic';
        [numnode1] = findnodeinnumddl(M,BC{i}.ddlbloque);
        [numnode2,ddl] = findnodeinnumddl(M,BC{i}.ddlfree);
        groupnode{i} = union(numnode1,numnode2);
    end
    ddl = unique(ddl);
    for k=1:length(ddl)
        ddl{k} = [ddl{k}, ','];
    end
    ddl = [ddl{:}];
    paramplot{i} = [ddl, name];
end

% [leg,I] = unique(paramplot);
% for i=1:length(BC)
%     [temp,loc] = ismember(paramplot{i},leg);
%     groupnodeplot{loc} = union(groupnodeplot{loc},groupnode{i});
% end
leg = paramplot;
groupnodeplot = groupnode;

Handles = [];
for i=1:length(groupnodeplot)
    nodep = getnode(M.node,groupnodeplot{i});
    marker = getpointstyles(i);
    markersize = getcharin('MarkerSize',varargin,8);
    Htemp = plot(nodep,marker,varargin{:});
    Handles = [Handles,Htemp];
end

dim = getindim(M);
node = M.node;
nodecoord = double(getcoord(node));
if ~(dim==3 && all(nodecoord(:,3)==nodecoord(1,3)))
    axis image
end
axis off
if ~isempty(Handles)
    legend(Handles,leg{:});
end

if nargout>=1
    varargout{1} = Handles;
end
if nargout>=2
    varargout{2} = leg;
end
