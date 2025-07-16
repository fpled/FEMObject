function plotV(u,scanm,varargin)
% function plotV(u,scanm,varargin)
% affichage des modes de la PCRADIALMATRIX
% scanm : ensembles des modes à tracer
% on utilise plot(V{i},varargin{:})
% varargin : arguments de la fonction plot
%
% plotV(u,scanm,'name',nom de la variable affichee)
% ce nom sera indicé par le numéro de la cellule correspondante
%
%
% plotV(u,scanm,'manutext',manutext)
% pour afficher manuellement les noms des modes
% manutext est une cellule
% manutext{1} = [x,y] position du texte dans la figure (corrdonnees absolues de la figure)
% manutext{2} = 'name'  nom du mode (il sera indice par le numero du mode)
% manutext{3} = arguments (cell array contenant les options de la fonction  text)
% appel de text(pos,name,arg{:})

if isempty(scanm) | nargin==1
    scanm = 1:u.m;
end

modename = getcharin('modename',varargin);
if ~isempty(modename)
    plottitle = 1;
else
    modename = inputname(1);
end

varargin=delcharin('modename',varargin);
plottitle = ischarin('plottitle',varargin);
if plottitle
  varargin=delonlycharin('plottitle',varargin);  
end

manutext = getcharin('manutext',varargin);
varargin=delonlycharin('manutext',varargin);

m = length(scanm);

nl=getcharin('nl',varargin,floor(sqrt(m)));
nc=getcharin('nc',varargin,ceil(m/nl));
for i=1:m
    if ischarin('fact',varargin);
        fullsubplot(nl,nc,i,getcharin('fact',varargin));
    else
fullsubplot(nl,nc,i);
    end
V = u.V{scanm(i)};

plot(V,varargin{:});
axis off
if ~isempty(modename) & plottitle
    title([ modename '_{' num2str(scanm(i)) '}'],'FontSize',16)
elseif ~isempty(manutext)
    modename=manutext{2};
    pos = manutext{1};
    name = manutext{2};
   
    arg = manutext(3:end);
   
    text(pos{:},[ modename '_{' num2str(scanm(i)) '}'],arg{:});
end

end
