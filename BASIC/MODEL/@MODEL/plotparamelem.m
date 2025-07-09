function varargout = plotparamelem(M,choix,varargin)
% function [Handles,leg] = plotparamelem(M,choix,varargin)
% choix : parametres a afficher
% Handles : handles vers les patchs d'elements
% leg : cellules contenant les legendes de chaque handle

node = getnode(M);
paramplot = {};
Handles = [];
dim = getindim(M);
optionsplot = patchoptions(dim,varargin{:});
leg = {};
for p=1:M.nbgroupelem
    elem = getgroupelem(M,p);
    switch choix
        case 'material'
            paramelem = getmaterialnumber(elem);
        case 'group'
            paramelem = p;
        otherwise
            paramelem = trygetparamelem(elem,choix);
    end
    
    if isempty(paramelem)
        optionsplot = setcharin('facecolor',optionsplot,'none');
        plot(elem,node,optionsplot{:});
        isalreadyplot = 1;
    else
        isalreadyplot = 0;
        for i=1:length(paramplot)
            if paramcmp(paramelem,paramplot{i})
                optionsplot = setcharin('facecolor',optionsplot,getfacecolor(i));
                plot(elem,node,optionsplot{:});
                isalreadyplot = 1;
                break
            end
        end
    end
    
    if isalreadyplot==0
        paramplot = [paramplot , {paramelem}];
        i = length(paramplot);
        
        if isa(paramelem,'double') || isa(paramelem,'logical')
            leg = [leg , {num2str(paramelem)}];
        elseif isa(paramelem,'char')
            leg = [leg , {paramelem}];
        else
            error('on ne peut pas afficher')
        end
        optionsplot = setcharin('facecolor',optionsplot,getfacecolor(i));
        Htemp = plot(elem,node,optionsplot{:});
        Handles = [Handles,Htemp];
    end
end
axis off
axis image

numview = getcharin('view',varargin);
up_vector = getcharin('camup',varargin);
camera_position = getcharin('campos',varargin);
if ~isempty(numview)
    view(numview)
elseif dim==3
    view(3)
end
if ~isempty(up_vector)
    camup(up_vector)
end
if ~isempty(camera_position)
    campos(camera_position);
end

if ~isempty(Handles)
    legend(Handles,leg{:});
end

if nargout>=1
    varargout{1} = Handles;
end
if nargout>=2
    varargout{2} = leg;
end


function rep =  paramcmp(a,b)
if isa(a,'logical')
    a = double(a);
end
if isa(b,'logical')
    b = double(b);
end

rep = strcmp(class(a),class(b));
if rep && isa(a,'char')
    rep = strcmp(a,b);
elseif rep && isa(a,'double')
    rep = all(size(a)==size(b));
    if rep
        rep = all(all(a==b));
    end
end
return

function paramelem = trygetparamelem(elem,choix)

try
    paramelem = getelementfield(elem,choix);
catch
    try
        paramelem = getelementgeomfield(elem,choix);
    catch
        try
            paramelem = getparam(elem,choix);
        catch
            paramelem = [];
        end
    end
end
return

