function plottext(u,text,varargin)
% function plottext(u,text)
% text : vecteur de caracteres ou de double
%    contenant une info a afficher pour chaque noeud
%
% function plottext(u,text,'Color',color,'FontSize',fontsize)

color = getcharin('Color',varargin,'b');
fontsize = getcharin('FontSize',varargin,14); 
plottext(u.POINT,text,'FontSize',fontsize,'Color',color);

