function optionsout=patchoptions(dim,varargin)
% function options=patchoptions(dim,indim,'propertyname',propertyvalue,...)
% dim : dimension
% propertyname       propertyvalue                          effet
% 'facecolor'        'k', ..., 'flat', 'interp'        voir patch
% 'edgecolor'        'k', ..., 'flat', 'interp'        voir patch
% 'linewidth'        a                      taille des edges (ne fonctionne pas avec edgecolor interp) 

if nargin==2 && isa(varargin{1},'double')
    indim=varargin{1};
else
    indim=dim;
end

facecolor = getcharin('facecolor',varargin,'none');
edgecolor = getcharin('edgecolor',varargin,'k');
optionsout = {'facecolor',facecolor,'edgecolor',edgecolor};
facevertexcdata = getcharin('facevertexcdata',varargin);
if indim==3 && ~ischarin('solid',varargin)
    % facealpha = getcharin('facealpha',varargin,0.3);
    facealpha = getcharin('facealpha',varargin,1);
    facelighting = getcharin('facelighting',varargin,'gouraud');
else
    facealpha = getcharin('facealpha',varargin);
    facelighting = getcharin('facelighting',varargin);
end
if ischarin('surface',varargin)
    optionsout = [optionsout, {'surface'}];
end
if ischarin('surfacemesh',varargin)
    optionsout = [optionsout, {'surfacemesh'}];
end

edgealpha = getcharin('edgealpha',varargin);
edgelighting = getcharin('edgelighting',varargin);    

if ~isempty(facevertexcdata)
    optionsout = setcharin('facevertexcdata', optionsout, facevertexcdata);    
end
if ~isempty(facealpha)
    optionsout = setcharin('facealpha', optionsout, facealpha);    
end
if ~isempty(edgealpha)
    optionsout = setcharin('edgealpha', optionsout, edgealpha);    
end
if ~isempty(facelighting)
    optionsout = setcharin('facelighting', optionsout, facelighting);    
end
if ~isempty(edgelighting)
    optionsout = setcharin('edgelighting', optionsout, edgelighting);    
end

if ischarin('noedges',varargin) && ~ischarin('edgecolor',varargin)
    optionsout = setcharin('edgecolor', optionsout, 'none');   
end
linewidth = getcharin('linewidth',varargin);
if linewidth
    optionsout = setcharin('linewidth', optionsout, linewidth);
end
marker = getcharin('marker',varargin);
if isempty(marker) && dim==0
    marker = '.';
end
if ~isempty(marker)
    optionsout = setcharin('marker', optionsout, marker);    
end
if ischarin('node',varargin)
    optionsout = setcharin('marker', optionsout, '.');    
end
markersize = getcharin('markersize',varargin);
if isempty(markersize) && dim==0
    markersize = 10;
end
if ~isempty(markersize)
    optionsout = setcharin('markersize', optionsout, markersize);    
end
markerfacecolor = getcharin('markerfacecolor',varargin);
markeredgecolor = getcharin('markeredgecolor',varargin);
if isempty(markerfacecolor) && dim==0
    markerfacecolor = edgecolor;
end
if isempty(markeredgecolor) && dim==0
    markeredgecolor = edgecolor;
end

if ~isempty(markerfacecolor)
    optionsout = setcharin('markerfacecolor', optionsout, markerfacecolor);    
end
if ~isempty(markeredgecolor)
    optionsout = setcharin('markeredgecolor', optionsout, markeredgecolor);    
end
