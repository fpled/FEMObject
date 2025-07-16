function varargout = lsenrichplot(M,varargin)

Handles = [];
leg = {};

plot(M)
num = getnumgroupelemwithfield(M,'lsenrich',1);
Htemp = plot(M,'selgroup',num,'FaceColor','y');
if length(Htemp)>0
legtemp = cell(1,length(Htemp(1)));
legtemp(:)={'enriched'};
Handles = [Handles,Htemp(1)];
leg  = [leg,legtemp];
end

   node=getnode(M);
   lsenrichnode = getlsenrichnode(node);
   repnature = getrepnature(node);
   enrichnature = getlsenrichnature(node);
   
   for j=1:length(enrichnature)
   num = find(repnature==j);
   nodep = getnode(node,lsenrichnode(num));
   Htemp = plot(nodep,getpointstyles(j));
   if length(Htemp)>0
   if ~isa(enrichnature{j},'char')
       enrichnature{j}=num2str(enrichnature{j});
   end
   leg = [leg,{enrichnature{j}}];
   Handles = [Handles,Htemp];
   end
   
   end
   
   legend(Handles,leg{:});

if nargout>=1
    varargout{1} = Handles;
end
if nargout==2
    varargout{2} = leg;
end
