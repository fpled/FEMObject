function varargout = gmshDomainWithSingleEdgeCrack(D,C,clD,clC,filename,indim,varargin)
% function varargout = gmshDomainWithSingleEdgeCrack(D,C,clD,clC,filename,indim)
% D : DOMAIN
% C : LINE in dim 2, QUADRANGLE in dim 3
% clD, clC : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getindim(D) by default)

Box         = getcharin('Box',varargin,[]);
noduplicate = ischarin('noduplicate',varargin);
refinecrack = ischarin('refinecrack',varargin);
recombine   = ischarin('recombine',varargin);

varargin = delcharin('Box',varargin);
varargin = delonlycharin({'noduplicate','refinecrack','recombine'},varargin);

if nargin<6 || isempty(indim)
    indim = getindim(D);
end
if nargin<4 || isempty(clC)
    clC = clD;
end

dim = getdim(D);

PD = getvertices(D);

if indim==2
    PC = getvertices(C);
    C = LINE(min(PC{:}),max(PC{:}));
    if refinecrack
        G = gmshfile(C,clC,[2,1],1);
    else
        G = gmshfile(C,[clD clC],[2,1],1);
    end
    G = createpoints(G,PD,clD,3:6);
    G = createcontour(G,2:6,2:6,1);
    G = createplanesurface(G,1,1);
    G = embedcurveinsurface(G,1,1);
    
    if recombine
        G = recombinesurface(G,1);
    end
    
    if ~noduplicate
        physicalgroup = 1;
        openboundaryphysicalgroup = 1;
        G = createphysicalpoint(G,2,openboundaryphysicalgroup);
        G = createphysicalcurve(G,1,physicalgroup);
    end

    G = createphysicalsurface(G,1,1);
    
elseif indim==3
    if refinecrack
        G = gmshfile(C,clC,1:4,1:4,1,1);
    else
        G = gmshfile(C,[clD clC clC clD],1:4,1:4,1,1);
    end
    G = createpoints(G,PD,clD,5:12);
    G = createcontour(G,[1 8 7 6 5],5:9,2);
    G = createplanesurface(G,2,2);
    G = embedcurveinsurface(G,1,2);
    
    G = createcontour(G,[4 9 10 11 12],10:14,3);
    G = createplanesurface(G,3,3);
    G = embedcurveinsurface(G,3,3);
    
    G = createlines(G,[[7 11];[12 8]],15:16);
    G = createcurveloop(G,-[6 16 13 15],4);
    G = createplanesurface(G,4,4);
    
    G = createlines(G,[[6 10];[9 5]],17:18);
    G = createcurveloop(G,[-8 17 -11 18],5);
    G = createplanesurface(G,5,5);
    
    G = createcurveloop(G,[-9 -18 -10 -14 16 -5],6);
    G = createplanesurface(G,6,6);
    G = embedcurveinsurface(G,4,6);
    
    G = createcurveloop(G,[15 -12 -17 -7],7);
    G = createplanesurface(G,7,7);
    
    if recombine
        G = recombinesurface(G);
    end
    
    G = createsurfaceloop(G,2:7,1);
    G = createvolume(G,1,1);
    G = embedsurfaceinvolume(G,1,1);
    
    if ~noduplicate
        physicalgroup = 1;
        openboundaryphysicalgroup = 1;
        G = createphysicalpoint(G,[1,4],openboundaryphysicalgroup);
        G = createphysicalcurve(G,[1,3,4],openboundaryphysicalgroup);
        G = createphysicalsurface(G,1,physicalgroup);
    end
    
    G = createphysicalvolume(G,1,1);
    
end

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

if nargin>=5 && ischar(filename)
    G = setfile(G,filename);
end

n = max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});

if ~noduplicate
    G = createcrack(G,dim-1,physicalgroup,openboundaryphysicalgroup);
    G = remesh(G,dim,varargin{:});
    
    [varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});
end
