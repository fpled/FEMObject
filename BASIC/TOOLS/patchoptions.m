function optionsout=patchoptions(dim,varargin)
% function options=patchoptions(dim,indim,'propertyname',propertyvalue,...)
% dim : dimension
% propertyname       propertyvalue                          effet
% 'FaceColor'        'k', ..., 'flat', 'interp'        voir patch
% 'EdgeColor'        'k', ..., 'flat', 'interp'        voir patch
% 'LineWidth'        a                      taille des edges (ne fonctionne pas avec edgecolor interp) 

if nargin==2 && isa(varargin{1},'double')
    indim=varargin{1};
else
    indim=dim;
end

facecolor = getcharin('FaceColor',varargin,'none');
edgecolor = getcharin('EdgeColor',varargin,'k');
optionsout = {'FaceColor',facecolor,'EdgeColor',edgecolor};
facevertexcdata = getcharin('FaceVertexCData',varargin);
if indim==3 && ~ischarin('solid',varargin)
    % facealpha = getcharin('FaceAlpha',varargin,0.3);
    facealpha = getcharin('FaceAlpha',varargin,1);
    facelighting = getcharin('FaceLighting',varargin,'gouraud');
else
    facealpha = getcharin('FaceAlpha',varargin);
    facelighting = getcharin('FaceLighting',varargin);
end
if ischarin('surface',varargin)
    optionsout = [optionsout, {'surface'}];
end
if ischarin('surfacemesh',varargin)
    optionsout = [optionsout, {'surfacemesh'}];
end

edgealpha = getcharin('EdgeAlpha',varargin);
edgelighting = getcharin('EdgeLighting',varargin);    

if ~isempty(facevertexcdata)
    optionsout = setcharin('FaceVertexCData', optionsout, facevertexcdata);    
end
if ~isempty(facealpha)
    optionsout = setcharin('FaceAlpha', optionsout, facealpha);    
end
if ~isempty(edgealpha)
    optionsout = setcharin('EdgeAlpha', optionsout, edgealpha);    
end
if ~isempty(facelighting)
    optionsout = setcharin('FaceLighting', optionsout, facelighting);    
end
if ~isempty(edgelighting)
    optionsout = setcharin('EdgeLighting', optionsout, edgelighting);    
end

if ischarin('noedges',varargin) && ~ischarin('EdgeColor',varargin)
    optionsout = setcharin('EdgeColor', optionsout, 'none');   
end
linewidth = getcharin('LineWidth',varargin);
if linewidth
    optionsout = setcharin('LineWidth', optionsout, linewidth);
end
marker = getcharin('Marker',varargin);
if isempty(marker) && dim==0
    marker = '.';
end
if ~isempty(marker)
    optionsout = setcharin('Marker', optionsout, marker);    
end
if ischarin('node',varargin)
    optionsout = setcharin('Marker', optionsout, '.');    
end
markersize = getcharin('MarkerSize',varargin);
if isempty(markersize) && dim==0
    markersize = 10;
end
if ~isempty(markersize)
    optionsout = setcharin('MarkerSize', optionsout, markersize);    
end
markerfacecolor = getcharin('MarkerFaceColor',varargin);
markeredgecolor = getcharin('MarkerEdgeColor',varargin);
if isempty(markerfacecolor) && dim==0
    markerfacecolor = edgecolor;
end
if isempty(markeredgecolor) && dim==0
    markeredgecolor = edgecolor;
end

if ~isempty(markerfacecolor)
    optionsout = setcharin('MarkerFaceColor', optionsout, markerfacecolor);    
end
if ~isempty(markeredgecolor)
    optionsout = setcharin('MarkerEdgeColor', optionsout, markeredgecolor);    
end
