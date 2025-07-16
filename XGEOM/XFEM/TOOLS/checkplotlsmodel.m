function checkplotlsmodel(S,varargin)
if isenrich(S)
    nl=3;
    nc=3;
else
    nl=2;
    nc=2;
end

set(gca,'ActivePositionProperty','OuterPosition')
clf
i=0;

i=i+1;subplot(nl,nc,i);
plotparamelem(S,'group')
title('groupelem')

i=i+1;subplot(nl,nc,i)
lstypeplot(S);
contourplot(getlevelsets(S),S,'LineWidth',2,'Color','m');
title('lstype')

i=i+1;subplot(nl,nc,i)
plotparamelem(S,'lsnumber')
title('lsnumber')

i=i+1;subplot(nl,nc,i)
plotparamelem(S,'lsnature')
title('lsnature')

if isenrich(S)
    i=i+1;subplot(nl,nc,i)
    plotparamelem(S,'lsenrich')
    title('lsenrich')
    
    i=i+1;subplot(nl,nc,i)
    lsenrichplot(S)
    contourplot(getlevelsets(S),S,'LineWidth',2,'Color','m');
    title('enrichment')
    
    i=i+1;subplot(nl,nc,i)
    plotparamnode(S,'nbddl','MarkerSize',8)
    title('nbddlpernode')
    
    if islscrackin(getlevelsets(S))
        i=i+1;subplot(nl,nc,i)
        plotparamelem(S,'tipnumber')
        title('tipnumber')
    end
    
end
