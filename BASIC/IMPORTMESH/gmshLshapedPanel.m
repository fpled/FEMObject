function varargout = gmshLshapedPanel(a,b,t,clD,clC,filename,indim,varargin)
% function varargout = gmshLshapedPanel(a,b,t,clD,clC,filename,indim)
% a : half-length
% b : distance of applied load from the right edge
% t : thickness
% clD, clC : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, 2 by default)

Box = getcharin('Box',varargin,[]);
varargin = delcharin('Box',varargin);

if nargin<7 || isempty(indim)
    indim = 2;
end
if nargin<5 || isempty(clC)
    clC = clD;
end

dim = indim;

% if indim==2
%     P{1} = [  0,  0];
%     P{2} = [a-b,  0];
%     P{3} = [  a,  0];
%     P{4} = [  a,  a];
%     P{5} = [ -a,  a];
%     P{6} = [ -a, -a];
%     P{7} = [  0, -a];
% elseif indim==3
    P{1} = [  0,  0, 0];
    P{2} = [a-b,  0, 0];
    P{3} = [  a,  0, 0];
    P{4} = [  a,  a, 0];
    P{5} = [ -a,  a, 0];
    P{6} = [ -a, -a, 0];
    P{7} = [  0, -a, 0];
% end

br = @(tag,k) sprintf('%s[%d]', tag, k); % bracket reference helper

G = GMSHFILE();
if nargin>=2 && ischar(filename)
    G = setfile(G,filename);
end

numpoints = 1:7;
numlines = 1:7;
numlineloop = 1;
numsurface = 1;
G = createpoints(G,P(1),clC,numpoints(1));
G = createpoints(G,P(2:7),clD,numpoints(2:7));
% G = createpoints(G,P(1:3),clC,numpoints(1:3));
% G = createpoints(G,P(4:7),clD,numpoints(4:7));
G = createcontour(G,numpoints,numlines,numlineloop);
G = createplanesurface(G,numlineloop,1);
if indim==3
    vect = [0,0,t];
    if ischarin('recombine',varargin) && ~ischarin('Layers',varargin)
        numlayers = max(1, round(t/min(clD,clC)));
        varargin = [varargin, {'Layers',numlayers}];
    end
    [G,tag] = extrude(G,vect,'Surface',numsurface,varargin{:});
    numtopsurface = br(tag,0); % top surface
    numvolume     = br(tag,1); % volume
end
if ischarin('recombine',varargin)
    if indim==2
        G = recombinesurface(G,numsurface);
    elseif indim==3
        G = recombinesurface(G,{numsurface,numtopsurface});
    end
end
if indim==2
    numphysicalsurface = 1;
    G = createphysicalsurface(G,numsurface,numphysicalsurface);
elseif indim==3
    numphysicalvolume = 1;
    G = createphysicalvolume(G,numvolume,numphysicalvolume);
end

varargin = delonlycharin('recombine',varargin);

% Box field
if ~isempty(Box) && isstruct(Box)
    if isfield(Box,'VIn')
        VIn = Box.VIn;
    else
        VIn = clC;
    end
    if isfield(Box,'VOut')
        VOut = Box.VOut;
    else
        VOut = clD;
    end
    XMin = Box.XMin;
    XMax = Box.XMax;
    YMin = Box.YMin;
    YMax = Box.YMax;
    if indim==3 || isfield(Box,'ZMin')
        ZMin = Box.ZMin;
    else
        ZMin = 0;
    end
    if indim==3 || isfield(Box,'ZMax')
        ZMax = Box.ZMax;
    else
        ZMax = 0;
    end
    if isfield(Box,'Thickness')
        Thickness = Box.Thickness;
    else
        Thickness = 0;
    end
    G = createboxfield(G,VIn,VOut,XMin,XMax,YMin,YMax,ZMin,ZMax,Thickness);
    G = setbgfield(G);
end

n = max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});
