function plotelemtype(ls,varargin)
% function plotelemtype(ls,varargin)

ls = actualise(ls,varargin{:});

groupout=[];
groupin=[];
groupcut=[];
for p=1:ls.D.nbgroupelem
    if all(lsisout(ls.D.groupelem{p},ls));
        groupout =[groupout,p];
    elseif all(lsisin(ls.D.groupelem{p},ls));
        groupin =[groupin,p];
    else
        groupcut =[groupcut,p];
    end
end

plot(ls.D,'selgroup',groupout,'FaceColor','r');
plot(ls.D,'group','selgroup',groupout,'FaceColor','g');
plot(ls.D,'group','selgroup',groupin,'FaceColor','y');
plot(ls.D,'group','selgroup',groupcut,'FaceColor','r');






