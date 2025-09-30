function varargout = gmshfilewithpoints_occ(D,P,clD,clP,numberembeddedpoints,number,varargin)
% function G = gmshfilewithpoints_occ(D,P,clD,clP,numberembeddedpoints,number)
% D : DOMAIN (rectangle in 2D or box in 3D)
% P : POINT
% clD, clP : characteristic length

extrusion = ischarin('extrude',varargin);
varargin = delonlycharin('extrude',varargin);

if ~iscell(P), P = {P}; end
if nargin<4 || isempty(clP), clP = clD; end
if isscalar(clP), clP = repmat(clP,1,length(P)); end

if nargin<5 || isempty(numberembeddedpoints)
    switch D.dim
        case 2
            numberembeddedpoints = 4+(1:length(P));
        case 3
            numberembeddedpoints = 8+(1:length(P));
    end
end
if nargin<6 || isempty(number), number = 1; end

P1 = getvertex(D,1);
switch D.dim
    case 2
        point = [P1,0];
    case 3
        point = P1;
end
extents = getsize(D);

br = @(tag,k) sprintf('%s[%d]', tag, k); % bracket reference helper

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
switch D.dim
    case 2
        G = createrectangle(G,point,extents,number);
    case 3
        if extrusion
            % Build base profile as a rectangle
            numbersurface = number;
            G = createrectangle(G,point,extents(1:2),numbersurface);
            % Extrude along the third direction using a translation vector
            % P5 = getvertex(D,5);
            % vect = P5-P1;
            vect = [0,0,extents(3)];
            if ischarin('recombine',varargin) && ~ischarin('Layers',varargin)
                numberlayers = max(1, round(extents(3)/clD));
                varargin = [varargin, {'Layers',numberlayers}];
            end
            [G,tag] = extrude(G,vect,'Surface',numbersurface,varargin{:});
            numbertopsurface  = br(tag,0); % top surface
            number            = br(tag,1); % volume
            numberlatsurfaces = arrayfun(@(k) br(tag,k), 1+(1:4), 'UniformOutput', false); % lateral surfaces
            numbersurfaces    = [{numbersurface, numbertopsurface}, numberlatsurfaces];
            if ischarin('recombine',varargin)
                G = recombinesurface(G,{numbersurface,numbertopsurface});
            end
        else
            G = createbox(G,point,extents,number);
        end
end
G = setmeshsize(G,clD);

G = createpoints(G,P,clP,numberembeddedpoints);
switch D.dim
    case 2
        G = embedpointsinsurface(G,numberembeddedpoints,number);
    case 3
        G = embedpointsinvolume(G,numberembeddedpoints,number);
end

if ischarin('recombine',varargin)
    if D.dim==2
        G = recombinesurface(G,number);
    elseif D.dim==3 && ~extrusion
        G = recombinesurface(G);
    end
end

varargout{1} = G;
varargout{2} = number;
if D.dim==3 && extrusion
    varargout{3} = numbersurfaces;
end
