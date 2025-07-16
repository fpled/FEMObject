function plottext(elem,node,text,varargin)

if (isa(text,'double') | isa(text,'char')) && length(text)==1
    xplot = calc_midpoint(elem,node);
else
    xplot = calc_midpointelem(elem,node);
end

color = getcharin('Color',varargin,'b');
fontsize = getcharin('FontSize',varargin,14); 

plottext(xplot,text,'Color',color,'FontSize',fontsize); 
