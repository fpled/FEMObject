function P =  PLOTOPTIONS(varargin)
% function P =  PLOTOPTIONS()
% contient les parametres pour les affichages graphiques
% permet d'executer des commandes graphiques sans les taper
% 

options = PARAMETERS('step',1,'clf',true,...
    'axison',false,'boxon',false,'boxstylefull',false,'axis',[],'caxis',[],'axisimage',false,'axissquare',false,...
    'view',[],'camup','auto','campos','auto','pause',true,...
    'colorbar',false,'colormap','default','fontsize',16,'makepause',false,'pausetime',0.01,...
    'xlim',[],'ylim',[],'zlim',[],'noxtick',false,'noytick',false,'noztick',false,'plotstep',1,'pbaspect',[],'daspect',[]); 

options = setparam(options,varargin{:});
P = struct();
P = class(P,'PLOTOPTIONS',options);
