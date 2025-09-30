function varargout = gmshfile(D,cl,numberpoints,numberlines,numberlineloop,numbersurface,varargin)
% function G = gmshfile(D,cl,numberpoints,numberlines,numberlineloop,numbersurface,numbersurfaceloop,numbervolume)
% D : DOMAIN (rectangle in 2D or box in 3D)
% cl : characteristic length

extrusion = ischarin('extrude',varargin);
varargin = delonlycharin('extrude',varargin);

switch D.dim
    case 2
        if nargin<3 || isempty(numberpoints), numberpoints = 1:4; end
        if nargin<4 || isempty(numberlines), numberlines = 1:4; end
        if nargin<5 || isempty(numberlineloop), numberlineloop = 1; end
        if nargin<6, numbersurface = []; elseif isempty(numbersurface), numbersurface = 1; end
    case 3
        if extrusion
            if nargin<3 || isempty(numberpoints), numberpoints = 1:4; end
            if nargin<4 || isempty(numberlines), numberlines = 1:4; end
            if nargin<5 || isempty(numberlineloop), numberlineloop = 1; end
            if nargin<6 || isempty(numbersurface), numbersurface = 1; end
        else
            if nargin<3 || isempty(numberpoints), numberpoints = 1:8; end
            if nargin<4 || isempty(numberlines), numberlines = 1:12; end
            if nargin<5 || isempty(numberlineloop), numberlineloop = 1:6; end
            if nargin<6 || isempty(numbersurface), numbersurface = 1:6; end
        end
        if nargin<7 || isempty(varargin) || ischar(varargin{1}), numbersurfaceloop = 1; else, numbersurfaceloop = varargin{1}; end
        if isempty(numbersurfaceloop), numbersurfaceloop = 1; end
        if nargin<7 || isempty(varargin) || ischar(varargin{1}) || (numel(varargin)>=2 && ischar(varargin{2})), numbervolume = 1; elseif nargin<8, numbervolume = []; else, numbervolume = varargin{2}; end
        % if isempty(numbervolume), numbervolume = 1; end
end

P = getvertices(D);

br = @(tag,k) sprintf('%s[%d]', tag, k); % bracket reference helper

G = GMSHFILE();
switch D.dim
    case 2
        G = createpoints(G,P,cl,numberpoints);
        G = createcontour(G,numberpoints,numberlines,numberlineloop);
        if ~isempty(numbersurface)
            G = createplanesurface(G,numberlineloop,numbersurface);
            if ischarin('recombine',varargin)
                G = recombinesurface(G,numbersurface);
            end
        end
    case 3
        if extrusion
            % Build base profile as a rectangle
            numberpoints = numberpoints(1:4);
            numberlines = numberlines(1:4);
            numberlineloop = numberlineloop(1);
            numbersurface = numbersurface(1);
            G = createpoints(G,P(1:4),cl,numberpoints);
            G = createcontour(G,numberpoints,numberlines,numberlineloop);
            G = createplanesurface(G,numberlineloop,numbersurface);
            % Extrude along the third direction using a translation vector
            vect = P{5}-P{1};
            if ischarin('recombine',varargin) && ~ischarin('Layers',varargin)
                numberlayers = max(1, round(vect(3)/cl));
                varargin = [varargin, {'Layers',numberlayers}];
            end
            [G,tag] = extrude(G,vect,'Surface',numbersurface,varargin{:});
            numbertopsurface  = br(tag,0); % top surface
            numbervolume      = br(tag,1); % volume
            numberlatsurfaces = arrayfun(@(k) br(tag,k), 1+(1:4), 'UniformOutput', false); % lateral surfaces
            numbersurfaces    = [{numbersurface, numbertopsurface}, numberlatsurfaces];
            if ischarin('recombine',varargin)
                G = recombinesurface(G,{numbersurface,numbertopsurface});
            end
        else
            G = createpoints(G,P(1:numel(numberpoints)),cl,numberpoints);
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
switch D.dim
    case 2
        varargout{2} = numbersurface;
    case 3
        varargout{2} = numbervolume;
        if extrusion
            varargout{3} = numbersurfaces;
        end
end
