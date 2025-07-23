function varargout = gmshfile(D,cl,numberpoints,numberlines,numberlineloop,numbersurface,varargin)
% function G = gmshfile(D,cl,numberpoints,numberlines,numberlineloop,numbersurface,numbersurfaceloop,numbervolume)
% D : DOMAIN (rectangle in 2D or box in 3D)
% cl : characteristic length

if D.dim==2
    if nargin<=2 || isempty(numberpoints), numberpoints = 1:4; end
    if nargin<=3 || isempty(numberlines), numberlines = 1:4; end
    if nargin<=4 || isempty(numberlineloop), numberlineloop = 1; numbersurface = 1; end
    if nargin==5, numbersurface = []; end
elseif D.dim==3
    if nargin<=2 || isempty(numberpoints), numberpoints = 1:8; end
    if nargin<=3 || isempty(numberlines), numberlines = 1:12; end
    if nargin<=4 || isempty(numberlineloop), numberlineloop = 1:6; end
    if nargin<=5 || isempty(numbersurface), numbersurface = 1:6; end
    if nargin<=6 || isempty(varargin{1}), numbersurfaceloop = 1; numbervolume = 1; end
    if nargin==7, numbersurfaceloop = varargin{1}; numbervolume = []; end
end

G = GMSHFILE();
P = getvertices(D);
if ischarin('extrude',varargin)
    G = createpoints(G,P(1:4),cl,numberpoints(1:4));
else
    G = createpoints(G,P,cl,numberpoints);
end
if D.dim==2
    G = createcontour(G,numberpoints,numberlines,numberlineloop);
    if ~isempty(numbersurface)
        G = createplanesurface(G,numberlineloop,numbersurface);
        if ischarin('recombine',varargin)
            G = recombinesurface(G,numbersurface);
        end
    end
elseif D.dim==3
    if ischarin('extrude',varargin)
        G = createcontour(G,numberpoints(1:4),numberlines(1:4),numberlineloop(1));
        G = createplanesurface(G,numberlineloop(1),numbersurface(1));
        vect = P{5}-P{1};
        [G,out] = extrude(G,vect,'Surface',numbersurface(1),varargin{:});
        numbertopsurface = [out,'[0]'];
        numbervolume     = [out,'[1]'];
        numberlateralsurfaces = {[out,'[2]'],[out,'[3]'],[out,'[4]'],[out,'[5]']};
        numbersurfaces = [{numbersurface}, {numbertopsurface}, numberlateralsurfaces{:}];
    else
        G = createdomaincontour(G,numberpoints,numberlines,numberlineloop,numbersurface,numbersurfaceloop);
        if ~isempty(numbervolume)
            G = createvolume(G,numbersurfaceloop,numbervolume);
        end
        if ischarin('recombine',varargin)
            G = recombinesurface(G);
        end
    end
end

varargout{1} = G;
if D.dim==2
    varargout{2} = numbersurface;
elseif D.dim==3
    varargout{2} = numbervolume;
    if ischarin('extrude',varargin)
        varargout{3} = numbersurfaces;
    end
end
