function multiplot(u,scanm,varargin)
% function multiplot(u,varargin)
% u : cell array contenant des entites � afficher
% 
% on utilise plot(u{i},varargin{:})
% varargin : arguments de la fonction plot
%
% multiplot(u,'name',nom de la variable affichee)
% ce nom sera indic� par le num�ro de la cellule correspondante
%
% multiplot(u,'manutext',manutext)
% pour afficher manuellement les noms des entites
% manutext est une cellule
% manutext{1} = [x,y] position du texte dans la figure (corrdonnees absolues de la figure)
% manutext{2} = 'name'  nom du mode (il sera indice par le numero du mode)
% manutext{3} = arguments (cell array contenant les options de la fonction  text)
% appel de text(pos,name,arg{:})



if isa(u,'double')
    warning('argument double, on interprete les colonnes comme des modes differents a afficher')    
    u = num2cell(u,1);
end
if ~isa(u,'cell')
    error('input argument must be a cell')
end

if nargin==1 || isempty(scanm)
    scanm = 1:length(u);
end

plottitle = ischarin('plottitle',varargin);
if plottitle
    varargin=delonlycharin('plottitle',varargin);  
end

modename = getcharin('name',varargin);
if ~isempty(modename)
    plottitle = 1;
else
    modename = inputname(1);
end

varargin=delcharin('name',varargin);

manutext = getcharin('manutext',varargin);
varargin=delonlycharin('manutext',varargin);

fulsub = ischarin('fullsubplot',varargin);
varargin = delonlycharin('fullsubplot',varargin);
m = length(scanm);
fact = getcharin('fact',varargin);
varargin = delcharin('fact',varargin);
nl=getcharin('nl',varargin,floor(sqrt(m)));
nc=getcharin('nc',varargin,ceil(m/nl));
varargin = delcharin('nl',varargin);
varargin = delcharin('nc',varargin);
for i=1:m
    if fulsub
        if ~isempty(fact)
            fullsubplot(nl,nc,i,fact);
        else
            fullsubplot(nl,nc,i);
        end
    else
        subplot(nl,nc,i);    
    end
    V = u{scanm(i)};

    plot(V,varargin{:});
    hold on
    if ~isempty(modename) && plottitle
        title([ modename '_{' num2str(scanm(i)) '}'],'FontSize',16)
    elseif ~isempty(manutext)
        modename=manutext{2};
        pos = manutext{1};
        arg = manutext(3:end);

        text(pos{:},[ modename '_{' num2str(scanm(i)) '}'],arg{:});
    end

end
