function varargout = gmshfilewithpoints(T,P,clT,clP,numbercenter,numberpoints,numberembeddedpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop,numbervolume,varargin)
% function G = gmshfilewithpoints(T,P,clT,clP,numbercenter,numberpoints,numberembeddedpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop,numbervolume)
% T : TORUS
% P : POINT
% clT, clP : characteristic length

extrusion = ischarin('extrude',varargin);
varargin = delonlycharin('extrude',varargin);

tol = getfemobjectoptions('tolerancepoint');
angle = T.angle;

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

if extrusion && ~isempty(angle_val) && angle_val <= tol
    error('gmshfilewithpoints:degenerateAngle', 'Angle must be > 0 for extrusion.');
end

if ~iscell(P), P = {P}; end
if nargin<4 || isempty(clP), clP = clT; end
if isscalar(clP) || isempty(clP), clP = repmat(clP,1,length(P)); end

% Generate all vertices (minor circle profile points)
PT = getvertices(T);
n = numel(PT)/4; % number of minor-circle rings in non-extruded mode

if extrusion
    if nargin<5 || isempty(numbercenter), numbercenter = 1:2; end                     % 1 major center + 1 minor center
    if nargin<6 || isempty(numberpoints), numberpoints = max(numbercenter)+(1:4); end % 4 points for the minor circle
    % if nargin<7 || isempty(numberembeddedpoints), numberembeddedpoints = max([numbercenter,numberpoints])+(1:length(P)); end
    if nargin<8 || isempty(numbercurves), numbercurves = 1:4; end                     % 4 arcs for the minor circle
    if nargin<9 || isempty(numbercurveloops), numbercurveloops = 1; end               % 1 minor circle
else
    if nargin<5 || isempty(numbercenter), numbercenter = 1:(n+3); end                   % 3 major centers + n minor centers
    if nargin<6 || isempty(numberpoints), numberpoints = max(numbercenter)+(1:4*n); end % 4 points for each of the n minor circles
    if nargin<7 || isempty(numberembeddedpoints), numberembeddedpoints = max([numbercenter,numberpoints])+(1:length(P)); end
    if nargin<8 || isempty(numbercurves)
        if isfull, numbercurves = 1:(8*n);   % 4*n minor arcs + 4*n major arcs
        else,      numbercurves = 1:(8*n-4); % 4*n minor arcs + 4*(n-1) major arcs
        end
    end
    if nargin<9 || isempty(numbercurveloops)
        if isfull, numbercurveloops = 1:(4*n);   % 4*n major circle arcs
        else,      numbercurveloops = 1:(4*n-2); % 4*(n-1) major circle arcs + 2 end minor circles
        end
    end
end
if nargin<10 || isempty(numbersurfaces), numbersurfaces = numbercurveloops; end
if nargin<11 || isempty(numbersurfaceloop), numbersurfaceloop = 1; end
if nargin<12 || isempty(numbervolume), numbervolume = 1; end

center = [T.cx,T.cy,T.cz];
radius = T.r1;
centers_minor = getcminor(T);

br = @(tag,k) sprintf('%s[%d]', tag, k); % bracket reference helper

G = GMSHFILE();
if extrusion
    center_minor = centers_minor(1,:);
    numbercenter = numbercenter(1:2);
    numberpoints = numberpoints(1:4);
    numbercurveloop = numbercurveloops(1);
    numbersurface   = numbersurfaces(1);
    G = createpoint(G,center,clT,numbercenter(1));
    G = createpoint(G,center_minor,clT,numbercenter(2));
    G = createpoints(G,PT(1:4),clT,numberpoints);
    G = createcirclecontour(G,numbercenter(2),numberpoints,numbercurves,numbercurveloop);
    G = createplanesurface(G,numbercurveloop,numbersurface);
    normal = [T.nx,T.ny,T.nz];
    axis = normal / norm(normal);
    if isfull
        % Full torus: closed major contour
        if ischarin('recombine',varargin) && ~ischarin('Layers',varargin)
            numberlayers = max(1, round((radius*pi/2)/clT));
            varargin = [varargin, {'Layers',numberlayers}];
        end
        numbertopsurfaces = cell(1,4);
        numbervolumes     = cell(1,4);
        numberlatsurfaces = cell(1,4);
        numberbasesurface = numbersurface;
        for i=1:4
            tag = ['out' num2str(i)];
            G = revolve(G,axis,center,'Pi/2','Surface',numberbasesurface,tag,varargin{:});
            numbertopsurfaces{i} = br(tag,0); % top surface
            numbervolumes{i}     = br(tag,1); % volume
            numberlatsurfaces{i} = arrayfun(@(k) br(tag,k), 1+(1:4), 'UniformOutput', false); % lateral surfaces
            numberbasesurface = numbertopsurfaces{i}; % base surface for next revolution
        end
        if ~ischarin('recombine',varargin)
            G = deletevolume(G,numbervolumes);
            G = deletesurface(G,numbertopsurfaces);
        end
        numbersurfaces = [numberlatsurfaces{:}];
        if ~ischarin('recombine',varargin)
            G = createsurfaceloop(G,numbersurfaces,numbersurfaceloop);
            G = createvolume(G,numbersurfaceloop,numbervolume);
        else
            numbervolume = numbervolumes;
            % G = recombinesurface(G,[{numbersurface}, numbertopsurfaces]);
            G = recombinesurface(G,numbertopsurfaces); % numbertopsurfaces{end} == numbersurface
        end
    else
        % Partial torus: open major circle arc with two end minor circles
        % Compute how many full Pi/2 chunks (N) and the remainder
        if isempty(angle_val)
            N = 0; % unknown numeric -> do one revolve with the char expression
        else
            N = min(3, max(0, floor((angle_val + tol)/(pi/2)))); % snap to exact quarters within tol, at most 3 quarter-turns
        end
        numbertopsurfaces = cell(1,N); % top surfaces from the Pi/2 chunks
        numbervolumes     = cell(1,N); % volumes from the Pi/2 chunks
        numberlatsurfaces = cell(1,N); % lateral surfaces from the Pi/2 chunks
        numberbasesurface = numbersurface;
        % N quarter-turns of Pi/2 (use char 'Pi/2' to keep exact GEO)
        args = varargin;
        if ischarin('recombine',varargin) && ~ischarin('Layers',varargin)
            numberlayers = max(1, round((radius*pi/2)/clT));
            args = [args, {'Layers',numberlayers}];
        end
        for i=1:N
            tag = ['out' num2str(i)];
            G = revolve(G,axis,center,'Pi/2','Surface',numberbasesurface,tag,args{:});
            numbertopsurfaces{i} = br(tag,0); % top surface
            numbervolumes{i}     = br(tag,1); % volume
            numberlatsurfaces{i} = arrayfun(@(k) br(tag,k), 1+(1:4), 'UniformOutput', false); % lateral surfaces
            numberbasesurface = numbertopsurfaces{i}; % base surface for next revolve (the new top surface)
        end
        % Remaining angle (expr + numeric when available)
        if ischar(angle)
            ang_expr = sprintf('(%s)-%d*Pi/2', angle_expr, N); % keep char math for exact GEO
        else
            ang_expr = angle - N*(pi/2);
        end
        if ~isempty(angle_val)
            ang_val = angle_val - N*(pi/2); % numeric if parsable
        else
            ang_val = []; % unknown numeric
        end
        % If remainder is (numerically) zero, skip a final revolve
        if ~isempty(ang_val) && abs(ang_val) <= tol
            % Delete the intermediate top surfaces only, keep the last quarter-turn top surface
            if N > 1
                if ~ischarin('recombine',varargin)
                    G = deletesurface(G,numbertopsurfaces(1:N-1));
                end
            end
            % Delete the temporary volumes from the Pi/2 chunks
            if N > 0
                if ~ischarin('recombine',varargin)
                    G = deletevolume(G,numbervolumes);
                end
                numbertopsurface  = numbertopsurfaces{N};   % last top surface
                numberlatsurfaces = [numberlatsurfaces{:}]; % flatten lateral surfaces to a 1x(4*N) cell array of char
            else
                % Degenerate case (angle == 0): both surfaces are the initial disk, no lateral surfaces
                numbertopsurface  = numbersurface;
                numberlatsurfaces = {};
            end
        else
            % Do the final (nonzero and non-quarter) revolve
            if ischarin('recombine',varargin) && ~ischarin('Layers',varargin) && ~isempty(ang_val)
                numberlayers = max(1, round((radius*ang_val)/clT));
                varargin = [varargin, {'Layers',numberlayers}];
            end
            [G,tag] = revolve(G,axis,center,ang_expr,'Surface',numberbasesurface,varargin{:});
            numbertopsurface      = br(tag,0); % top surface
            numberlastvolume      = br(tag,1); % last volume
            numberlastlatsurfaces = arrayfun(@(k) br(tag,k), 1+(1:4), 'UniformOutput', false); % last lateral surfaces
            numbervolumes     = [numbervolumes, {numberlastvolume}];
            numberlatsurfaces = [numberlatsurfaces{:}, numberlastlatsurfaces];
            if ~ischarin('recombine',varargin)
                G = deletevolume(G, numbervolumes);
                % Delete the intermediate top surfaces, keep the final one
                if ~isempty(numbertopsurfaces)
                    G = deletesurface(G,numbertopsurfaces);
                end
            end
        end
        numbersurfaces = [{numbersurface, numbertopsurface}, numberlatsurfaces];
        if ~ischarin('recombine',varargin)
            G = createsurfaceloop(G,numbersurfaces,numbersurfaceloop);
            G = createvolume(G,numbersurfaceloop,numbervolume);
        else
            numbervolume = numbervolumes;
            % G = recombinesurface(G,[{numbersurface}, numbertopsurfaces]);
            G = recombinesurface(G,numbertopsurfaces);
        end
    end
else
    centers_major = getcmajor(T);
    G = createpoint(G,center,clT,numbercenter(1));
    G = createpoints(G,centers_major,clT,numbercenter(2:3));
    G = createpoints(G,centers_minor,clT,numbercenter(4:end));
    G = createpoints(G,PT,clT,numberpoints);
    if isfull
        % Full torus
        G = createtoruscontour(G,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop);
    else
        % Partial torus
        G = createtorusarccontour(G,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop);
    end
    G = createvolume(G,numbersurfaceloop,numbervolume);
end

if extrusion
    if nargin<7 || isempty(numberembeddedpoints)
        tagembeddedpoint = 'p';
        G = newp(G,tagembeddedpoint);
        numberembeddedpoints = [{tagembeddedpoint}, arrayfun(@(k) sprintf('%s+%d', tagembeddedpoint, k), 1:(length(P)-1), 'UniformOutput', false)];
    end
end
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
