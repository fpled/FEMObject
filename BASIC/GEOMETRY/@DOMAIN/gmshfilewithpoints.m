function varargout = gmshfilewithpoints(D,P,clD,clP,numberpoints,numberembeddedpoints,numberlines,numberlineloop,numbersurface,varargin)
% function G = gmshfilewithpoints(D,P,clD,clP,numberpoints,numberembeddedpoints,numberlines,numberlineloop,numbersurface,numbersurfaceloop,numbervolume)
% D : DOMAIN (rectangle in 2D or box in 3D)
% P : POINT
% clD, clP : characteristic length

extrusion = ischarin('extrude',varargin);
varargin = delonlycharin('extrude',varargin);

if ~iscell(P), P = {P}; end
if nargin<4 || isempty(clP), clP = clD; end
if isscalar(clP), clP = repmat(clP,1,length(P)); end

switch D.dim
    case 2
        if nargin<5 || isempty(numberpoints), numberpoints = 1:4; end
        if nargin<6 || isempty(numberembeddedpoints), numberembeddedpoints = max(numberpoints)+(1:length(P)); end
        if nargin<7 || isempty(numberlines), numberlines = 1:4; end
        if nargin<8 || isempty(numberlineloop), numberlineloop = 1; end
        if nargin<9 || isempty(numbersurface), numbersurface = 1; end
    case 3
        if extrusion
            if nargin<5 || isempty(numberpoints), numberpoints = 1:4; end
            % if nargin<6 || isempty(numberembeddedpoints), numberembeddedpoints = max(numberpoints)+(1:length(P)); end
            if nargin<7 || isempty(numberlines), numberlines = 1:4; end
            if nargin<8 || isempty(numberlineloop), numberlineloop = 1; end
            if nargin<9 || isempty(numbersurface), numbersurface = 1; end
        else
            if nargin<5 || isempty(numberpoints), numberpoints = 1:8; end
            if nargin<6 || isempty(numberembeddedpoints), numberembeddedpoints = max(numberpoints)+(1:length(P)); end
            if nargin<7 || isempty(numberlines), numberlines = 1:12; end
            if nargin<8 || isempty(numberlineloop), numberlineloop = 1:6; end
            if nargin<9 || isempty(numbersurface), numbersurface = 1:6; end
        end
        if nargin<10 || isempty(varargin) || ischar(varargin{1}), numbersurfaceloop = 1; else, numbersurfaceloop = varargin{1}; end
        if isempty(numbersurfaceloop), numbersurfaceloop = 1; end
        if nargin<11 || isempty(varargin) || ischar(varargin{1}) || (numel(varargin)>=2 && ischar(varargin{2})), numbervolume = 1; else, numbervolume = varargin{2}; end
        % if isempty(numbervolume), numbervolume = 1; end
end

PD = getvertices(D);

br = @(tag,k) sprintf('%s[%d]', tag, k); % bracket reference helper

G = GMSHFILE();
switch D.dim
    case 2
        G = createpoints(G,PD,clD,numberpoints);
        G = createcontour(G,numberpoints,numberlines,numberlineloop);
        G = createplanesurface(G,numberlineloop,numbersurface);
    case 3
        if extrusion
            % Build base profile as a rectangle
            numberpoints = numberpoints(1:4);
            numberlines = numberlines(1:4);
            numberlineloop = numberlineloop(1);
            numbersurface = numbersurface(1);
            G = createpoints(G,PD(1:4),clD,numberpoints);
            G = createcontour(G,numberpoints,numberlines,numberlineloop);
            G = createplanesurface(G,numberlineloop,numbersurface);
            % Extrude along the third direction using a translation vector
            vect = PD{5}-PD{1};
            if ischarin('recombine',varargin) && ~ischarin('Layers',varargin)
                numberlayers = max(1, round(vect(3)/clD));
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
            G = createpoints(G,PD(1:numel(numberpoints)),clD,numberpoints);
            G = createdomaincontour(G,numberpoints,numberlines,numberlineloop,numbersurface,numbersurfaceloop);
            G = createvolume(G,numbersurfaceloop,numbervolume);
        end
end

if D.dim==3 && extrusion
    if nargin<6 || isempty(numberembeddedpoints)
        tagembeddedpoint = 'p';
        G = newp(G,tagembeddedpoint);
        numberembeddedpoints = [{tagembeddedpoint}, arrayfun(@(k) sprintf('%s+%d', tagembeddedpoint, k), 1:(length(P)-1), 'UniformOutput', false)];
    end
end
G = createpoints(G,P,clP,numberembeddedpoints);
switch D.dim
    case 2
        G = embedpointsinsurface(G,numberembeddedpoints,numbersurface);
    case 3
        G = embedpointsinvolume(G,numberembeddedpoints,numbervolume);
end

if ischarin('recombine',varargin)
    if D.dim==2
        G = recombinesurface(G,numbersurface);
    elseif D.dim==3 && ~extrusion
        G = recombinesurface(G);
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
