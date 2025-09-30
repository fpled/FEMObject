function varargout = gmshfilewithpoints_occ(C,P,clC,clP,numberembeddedpoints,numbervolume,varargin)
% function G = gmshfilewithpoints_occ(C,P,clC,clP,numberembeddedpoints,numbervolume)
% C : CYLINDER
% P : POINT
% clC, clP : characteristic length

extrusion = ischarin('extrude',varargin);
varargin = delonlycharin('extrude',varargin);

tol = getfemobjectoptions('tolerancepoint');
angle = C.angle;

% Normalize and (try to) evaluate angle
if isstring(angle), angle = char(angle); end
if ischar(angle)
    angle_expr = angle;
    angle_val  = str2num(lower(angle));
else
    angle_expr = num2str(angle);
    angle_val  = angle;
end

% Detect full/partial revolution when numeric, or by normalized char if not evaluable
if ~isempty(angle_val)
    isfull = abs(angle_val - 2*pi) < tol;
else
    s = regexprep(lower(angle_expr),'\s+',''); % ignore case and remove spaces (strip whitespace)
    isfull = strcmp(s,'2*pi');
end

if ~iscell(P), P = {P}; end
if nargin<4 || isempty(clP), clP = clC; end
if isscalar(clP), clP = repmat(clP,1,length(P)); end

if nargin<5 || isempty(numberembeddedpoints)
    if isfull, numberembeddedpoints = 2+(1:length(P));
    else,      numberembeddedpoints = 6+(1:length(P));
    end
end
if nargin<6 || isempty(numbervolume), numbervolume = 1; end

center = [C.cx,C.cy,C.cz];
radius = C.r;
height = C.h;

n = [C.nx,C.ny,C.nz];
v = [C.vx,C.vy];
normal = n / norm(n);
axis   = height * normal;
[dir,ang] = calcrotation_direction_angle(C,v,n);

br = @(tag,k) sprintf('%s[%d]', tag, k); % bracket reference helper

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
if extrusion
    % Build base profile as a disk (full) or plane surface = circle arc + radial lines (partial)
    if isfull
        % Full cylinder base: disk (closed circle contour)
        % numbercurve = 'c';
        % numbercurveloop = 'cl';
        numbersurface = 's';
        % G = newc(G,numbercurve);
        % G = createcircle(G,center,radius,numbercurve);
        % G = newcl(G,numbercurveloop);
        % G = createcurveloop(G,numbercurve,numbercurveloop);
        % G = news(G,numbersurface);
        % G = createplanesurface(G,numbercurveloop,numbersurface);
        G = news(G,numbersurface);
        G = createdisk(G,center,radius,numbersurface);
    else
        % Partial cylinder base: open circle arc + two radial lines
        numbercurve = 'c';
        numberpoints = 'p';
        numbercenter = 'pc';
        numberlines = {'ls','le'};
        numbercurves = {numberlines{1},numbercurve,numberlines{2}};
        numbercurveloop = 'cl';
        numbersurface = 's';
        if ischar(angle)
            angles = {0,angle};
        else
            angles = [0,angle];
        end
        G = newc(G,numbercurve);
        G = createcircle(G,center,radius,numbercurve,angles); % open circle arc
        G = pointsofcurve(G,numbercurve,numberpoints);
        G = newp(G,numbercenter);
        G = createpoint(G,center,clC,numbercenter);
        G = newc(G,numberlines{1});
        G = createline(G,{numbercenter,br(numberpoints,0)},numberlines{1}); % first (start) radial line
        G = newc(G,numberlines{2});
        G = createline(G,{br(numberpoints,1),numbercenter},numberlines{2}); % last (end) radial line
        G = newcl(G,numbercurveloop);
        G = createcurveloop(G,numbercurves,numbercurveloop);
        G = news(G,numbersurface);
        G = createplanesurface(G,numbercurveloop,numbersurface);
    end
    if abs(ang) > eps
        G = rotate(G,dir,center,ang,'Surface',numbersurface);
    end
    if ischarin('recombine',varargin) && ~ischarin('Layers',varargin)
        numberlayers = max(1, round(height/clC));
        varargin = [varargin, {'Layers',numberlayers}];
    end
    [G,tag] = extrude(G,axis,'Surface',numbersurface,varargin{:});
    numbertopsurface = br(tag,0); % top surface
    numbervolume     = br(tag,1); % volume
    numberlatsurface = br(tag,2); % lateral surface
    numbersurfaces   = {numbersurface, numbertopsurface, numberlatsurface};
    if ischarin('recombine',varargin)
        G = recombinesurface(G,{numbersurface,numbertopsurface});
    end
else
    % Use an arbitrary tangent vector v
    % if isfull
    %     % Full cylinder
    %     G = createcylinder(G,center,axis,radius,numbervolume);
    % else
    %     % Partial cylinder
    %     G = createcylinder(G,center,axis,radius,numbervolume,angle);
    % end
    % Use the predefined tangent vector v
    z = [0,0,1];
    zaxis = height * z;
    if isfull
        % Full cylinder
        G = createcylinder(G,center,zaxis,radius,numbervolume);
    else
        % Partial cylinder
        G = createcylinder(G,center,zaxis,radius,numbervolume,angle);
    end
    if abs(ang) > eps
        G = rotate(G,dir,center,ang,'Volume',numbervolume);
    end
end
G = setmeshsize(G,clC);

G = createpoints(G,P,clP,numberembeddedpoints);
G = embedpointsinvolume(G,numberembeddedpoints,numbervolume);

if ischarin('recombine',varargin) && ~extrusion
    G = recombinesurface(G);
end

varargout{1} = G;
varargout{2} = numbervolume;
if extrusion
    varargout{3} = numbersurfaces;
end
