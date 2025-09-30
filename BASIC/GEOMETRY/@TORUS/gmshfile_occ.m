function varargout = gmshfile_occ(T,cl,numbervolume,varargin)
% function G = gmshfile_occ(T,cl,numbervolume)
% T : TORUS
% cl : characteristic length

extrusion = ischarin('extrude',varargin);
varargin = delonlycharin('extrude',varargin);

tol = getfemobjectoptions('tolerancepoint');

if nargin<3 || isempty(numbervolume), numbervolume = 1; end

center = [T.cx,T.cy,T.cz];
radii  = [T.r1,T.r2];
angle  = T.angle;
centers_minor = getcminor(T);

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

n = [T.nx,T.ny,T.nz];
v = [T.vx,T.vy];
axis = n / norm(n);
[dir,ang] = calcrotation_direction_angle(T,v,n);

br = @(tag,k) sprintf('%s[%d]', tag, k); % bracket reference helper

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
if extrusion
    center_minor = centers_minor(1,:);
    radius_minor = radii(2); % radius_minor = T.r2;
    dir_minor = [1,0,0];
    ang_minor = 'Pi/2';
    if isfull
        % Full torus: closed major contour
        if ~ischarin('recombine',varargin)
            %% Curve extrusion
            numbercurve = 'c';
            % numbercurveloop = 'cl';
            % numbersurface = 's';
            G = newc(G,numbercurve);
            G = createcircle(G,center_minor,radius_minor,numbercurve);
            % G = newcl(G,numbercurveloop);
            % G = createcurveloop(G,numbercurve,numbercurveloop);
            % G = news(G,numbersurface);
            % G = createplanesurface(G,numbercurveloop,numbersurface);
            G = rotate(G,dir_minor,center_minor,ang_minor,'Curve',numbercurve);
            if abs(ang) > eps
                G = rotate(G,dir,center,ang,'Curve',numbercurve);
            end
            tag1 = 'out1'; tag2 = 'out2';
            G = revolve(G,axis,center,'Pi','Curve',numbercurve,tag1,varargin{:});
            G = revolve(G,axis,center,'-Pi','Curve',numbercurve,tag2,varargin{:});
            % numbertopcurves = {br(tag1,0), br(tag2,0)};
            numbersurfaces  = {br(tag1,1), br(tag2,1)};
            % numberlatcurves = {br(tag1,2), br(tag2,2)};
            
            %% Surface extrusion
            % numbersurface = 's';
            % G = news(G,numbersurface);
            % G = createdisk(G,center_minor,radius_minor,numbersurface);
            % G = rotate(G,dir_minor,center_minor,ang_minor,'Surface',numbersurface);
            % if abs(ang) > eps
            %     G = rotate(G,dir,center,ang,'Surface',numbersurface);
            % end
            % tag1 = 'out1'; tag2 = 'out2';
            % G = revolve(G,axis,center,'Pi','Surface',numbersurface,tag1,varargin{:});
            % G = revolve(G,axis,center,'-Pi','Surface',numbersurface,tag2,varargin{:});
            % numbertopsurfaces = {br(tag1,0), br(tag2,0)};
            % numbervolumes     = {br(tag1,1), br(tag2,1)};
            % numbersurfaces    = {br(tag1,2), br(tag2,2)};
            % G = deletevolume(G,numbervolumes);
            % G = deletesurface(G,[{numbersurface}, numbertopsurfaces]);
            
            G = setcoherence(G); % remove all duplicate elementary entities
            numbersurfaceloop = 'sl';
            G = newsl(G,numbersurfaceloop);
            G = createsurfaceloop(G,numbersurfaces,numbersurfaceloop);
            G = createvolume(G,numbersurfaceloop,numbervolume);
            G = deletesurface(G,numbersurfaces);
        else
            %% Surface extrusion
            numbersurface = 's';
            G = news(G,numbersurface);
            G = createdisk(G,center_minor,radius_minor,numbersurface);
            G = rotate(G,dir_minor,center_minor,ang_minor,'Surface',numbersurface);
            if abs(ang) > eps
                G = rotate(G,dir,center,ang,'Surface',numbersurface);
            end
            if ischarin('recombine',varargin) && ~ischarin('Layers',varargin)
                numberlayers = max(1, round((radii(1)*pi)/cl));
                varargin = [varargin, {'Layers',numberlayers}];
            end
            tag1 = 'out1'; tag2 = 'out2';
            G = revolve(G,axis,center,'Pi','Surface',numbersurface,tag1,varargin{:});
            G = revolve(G,axis,center,'-Pi','Surface',numbersurface,tag2,varargin{:});
            % numbertopsurfaces = {br(tag1,0), br(tag2,0)};
            numbervolume      = {br(tag1,1), br(tag2,1)};
            numbersurfaces    = {br(tag1,2), br(tag2,2)};
            G = setcoherence(G); % remove all duplicate elementary entities
            G = recombinesurface(G);
        end
    else
        % Partial torus: open major circle arc with two end minor circles
        % numbercurve = 'c';
        % numbercurveloop = 'cl';
        % numbersurface = 's';
        % G = newc(G,numbercurve);
        % G = createcircle(G,center_minor,radius_minor,numbercurve);
        % G = newcl(G,numbercurveloop);
        % G = createcurveloop(G,numbercurve,numbercurveloop);
        % G = news(G,numbersurface);
        % G = createplanesurface(G,numbercurveloop,numbersurface);
        
        numbersurface = 's';
        G = news(G,numbersurface);
        G = createdisk(G,center_minor,radius_minor,numbersurface);
        
        G = rotate(G,dir_minor,center_minor,ang_minor,'Surface',numbersurface);
        if abs(ang) > eps
            G = rotate(G,dir,center,ang,'Surface',numbersurface);
        end
        if ischarin('recombine',varargin) && ~ischarin('Layers',varargin) && ~isempty(angle_val)
            numberlayers = max(1, round((radii(1)*angle_val)/cl));
            varargin = [varargin, {'Layers',numberlayers}];
        end
        [G,tag] = revolve(G,axis,center,angle,'Surface',numbersurface,varargin{:});
        numbertopsurface = br(tag,0);
        numbervolume     = br(tag,1);
        numberlatsurface = br(tag,2);
        numbersurfaces = {numbersurface, numbertopsurface, numberlatsurface};
        if ischarin('recombine',varargin)
            G = recombinesurface(G,{numbersurface,numbertopsurface});
        end
    end
else
    if isfull
        % Full torus
        G = createtorus(G,center,radii,numbervolume);
    else
        % Partial torus
        G = createtorus(G,center,radii,numbervolume,angle);
    end
    if abs(ang) > eps
        G = rotate(G,dir,center,ang,'Volume',numbervolume);
    end
    if ischarin('recombine',varargin)
        G = recombinesurface(G);
    end
end
G = setmeshsize(G,cl);

varargout{1} = G;
varargout{2} = numbervolume;
if extrusion
    varargout{3} = numbersurfaces;
end
