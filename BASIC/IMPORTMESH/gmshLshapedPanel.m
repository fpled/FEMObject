function varargout = gmshLshapedPanel(a,b,t,clD,clC,filename,indim,varargin)
% function varargout = gmshLshapedPanel(a,b,t,clD,clC,filename,indim)
% a : half-length
% b : distance of applied load from the right edge
% t : thickness
% clD, clC : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, 2 by default)

if nargin<7 || isempty(indim)
    indim = 2;
end
if nargin<5 || isempty(clC)
    clC = clD;
end

if indim==2
    P{1} = [  0,  0];
    P{2} = [a-b,  0];
    P{3} = [  a,  0];
    P{4} = [  a,  a];
    P{5} = [ -a,  a];
    P{6} = [ -a, -a];
    P{7} = [  0, -a];
elseif indim==3
    P{1} = [  0,  0, 0];
    P{2} = [a-b,  0, 0];
    P{3} = [  a,  0, 0];
    P{4} = [  a,  a, 0];
    P{5} = [ -a,  a, 0];
    P{6} = [ -a, -a, 0];
    P{7} = [  0, -a, 0];
end

G = GMSHFILE();
if nargin>=2 && ischar(filename)
    G = setfile(G,filename);
end

% G = createpoints(G,P(1),clC,1);
% G = createpoints(G,P(2:7),clD,2:7);
G = createpoints(G,P(1:3),clC,1:3);
G = createpoints(G,P(4:7),clD,4:7);
G = createcontour(G,1:7,1:7,8);
G = createplanesurface(G,8,1);
if indim==3
    vect = [0,0,t];
    [G,out] = extrude(G,vect,'Surface',1,varargin{:});
    numbervolume = [out,'[1]'];
end
if ischarin('recombine',varargin)
    if indim==2
        G = recombinesurface(G,1);
    elseif indim==3
        G = recombinesurface(G);
    end
end
if indim==2
    G = createphysicalsurface(G,1,1);
elseif indim==3
    G = createphysicalvolume(G,numbervolume,1);
end

varargin = delonlycharin('recombine',varargin);

% Box field
B = getcharin('Box',varargin,[]);
if ~isempty(B) && isstruct(B)
    if isfield(B,'VIn')
        VIn = B.VIn;
    else
        VIn = clC;
    end
    if isfield(B,'VOut')
        VOut = B.VOut;
    else
        VOut = clD;
    end
    XMin = B.XMin;
    XMax = B.XMax;
    YMin = B.YMin;
    YMax = B.YMax;
    if indim==3 || isfield(B,'ZMin')
        ZMin = B.ZMin;
    else
        ZMin = 0;
    end
    if indim==3 || isfield(B,'ZMax')
        ZMax = B.ZMax;
    else
        ZMax = 0;
    end
    if isfield(B,'Thickness')
        Thickness = B.Thickness;
    else
        Thickness = 0;
    end
    G = createboxfield(G,VIn,VOut,XMin,XMax,YMin,YMax,ZMin,ZMax,Thickness);
    G = setbgfield(G);
end

n = max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,indim:-1:indim-n+1,varargin{:});
