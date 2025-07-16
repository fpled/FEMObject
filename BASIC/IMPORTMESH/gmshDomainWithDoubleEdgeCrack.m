function varargout = gmshDomainWithDoubleEdgeCrack(D,Ca,Cb,clD,clC,filename,indim,varargin)
% function varargout = gmshDomainWithDoubleEdgeCrack(D,Ca,Cb,clD,clC,filename,indim)
% D : DOMAIN
% Ca, Cb : LIGNE in dim 2, QUADRANGLE in dim 3
% clD, clC : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getindim(D) by default)

noduplicate = ischarin('noduplicate',varargin);
varargin = delonlycharin('noduplicate',varargin);

if nargin<7 || isempty(indim)
    indim = getindim(D);
end
if nargin<5 || isempty(clC)
    clC = clD;
end

PD = getvertices(D);

if indim==2
    PCa = getvertices(Ca);
    PCb = getvertices(Cb);
    Ca = LIGNE(min(PCa{:}),max(PCa{:}));
    Cb = LIGNE(min(PCb{:}),max(PCb{:}));
    if ischarin('refinecrack',varargin)
        Ga = gmshfile(Ca,clC,[2 1],1);
        Gb = gmshfile(Cb,clC,[8 5],8);
    else
        Ga = gmshfile(Ca,[clD clC],[2 1],1);
        Gb = gmshfile(Cb,[clD clC],[8 5],8);
    end
    G = Ga+Gb;
    G = createpoints(G,PD,clD,[3:4,6:7]);
    G = createcontour(G,2:7,2:7,1);
    G = createplanesurface(G,1,1);
    G = embedcurvesinsurface(G,[1,8],1);
    if ~noduplicate
        physicalgroup = 1;
        openboundaryphysicalgroup = 1;
        G = createphysicalpoint(G,[2,5],openboundaryphysicalgroup);
        G = createphysicalcurve(G,[1,8],physicalgroup);
    end
    if ischarin('recombine',varargin)
        G = recombinesurface(G,1);
    end
    G = createphysicalsurface(G,1,1);
    
elseif indim==3
    if ischarin('refinecrack',varargin)
        Ga = gmshfile(Ca,clC,1:4,1:4,1,1);
        Gb = gmshfile(Cb,clC,5:8,5:8,2,2);
    else
        Ga = gmshfile(Ca,[clD clC clC clD],1:4,1:4,1,1);
        Gb = gmshfile(Cb,[clD clC clC clD],5:8,5:8,2,2);
    end
    G = Ga+Gb;
    G = createpoints(G,PD,clD,9:16);
    G = createcontour(G,[1 12 11 5 10 9],9:14,3);
    G = createplanesurface(G,3,3);
    G = embedcurvesinsurface(G,[1 5],3);
    
    G = createcontour(G,[4 13 14 8 15 16],15:20,4);
    G = createplanesurface(G,4,4);
    G = embedcurvesinsurface(G,[3 7],4);
    
    G = createlines(G,[[11 15];[16 12]],21:22);
    G = createcurveloop(G,-[10 22 19 21],5);
    G = createplanesurface(G,5,5);
    
    G = createlines(G,[[10 14];[13 9]],23:24);
    G = createcurveloop(G,[-13 23 -16 24],6);
    G = createplanesurface(G,6,6);
    
    G = createcurveloop(G,[-14 -24 -15 -20 22 -9],7);
    G = createplanesurface(G,7,7);
    G = embedcurveinsurface(G,4,7);
    
    G = createcurveloop(G,[-11 21 -18 -17 -23 -12],8);
    G = createplanesurface(G,8,8);
    G = embedcurveinsurface(G,8,8);
    
    if ischarin('recombine',varargin)
        G = recombinesurface(G,1);
        G = recombinesurface(G,2);
        G = recombinesurface(G,3);
        G = recombinesurface(G,4);
        G = recombinesurface(G,5);
        G = recombinesurface(G,6);
        G = recombinesurface(G,7);
    end
    G = createsurfaceloop(G,3:8,1);
    G = createvolume(G,1,1);
    G = embedsurfacesinvolume(G,[1 2],1);
    if ~noduplicate
        physicalgroup = 1;
        openboundaryphysicalgroup = 1;
        G = createphysicalpoint(G,[1 4 5 8],openboundaryphysicalgroup);
        G = createphysicalcurve(G,[1 3 4 5 6 7],openboundaryphysicalgroup);
        G = createphysicalsurface(G,1,physicalgroup);
    end
    G = createphysicalvolume(G,1,1);
    
end
varargin = delonlycharin({'recombine','refinecrack'},varargin);

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

if nargin>=6 && ischar(filename)
    G = setfile(G,filename);
end

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(D):-1:getdim(D)-n+1,varargin{:});

if ~noduplicate
    G = createcrack(G,getdim(D)-1,physicalgroup,openboundaryphysicalgroup);
    G = remesh(G,getdim(D),varargin{:});
    G = deleteoptfile(G);
    
    [varargout{:}] = gmsh2femobject(indim,G,getdim(D):-1:getdim(D)-n+1,varargin{:});
end
