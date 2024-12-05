function varargout = gmshdomainwithedgenotch(D,C,c,clD,clC,filename,indim,varargin)
% function varargout = gmshdomainwithedgenotch(D,C,c,clD,clC,filename,indim,notchtype)
% D : DOMAIN
% C : LIGNE in dim 2, QUADRANGLE in dim 3
% c : width of the edge notch
% clD, clC : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getindim(D) by default)
% notchtype : 'c', 'r' or 'v' for circular, rectangular or V notch (optional, 'c' by default)

if nargin<7 || isempty(indim)
    indim = getindim(D);
end
if nargin<5 || isempty(clC)
    clC = clD;
end

PD = getvertices(D);
PC = getvertices(C);

if indim==2
    PB{1} = PC{1} + [0,c/2];
    PB{2} = PC{1} - [0,c/2];
    G = GMSHFILE();
    if ischarin('refinecrack',varargin)
        clB = clC;
    else
        clB = clD;
    end
    G = createpoints(G,PB,clB,[7,2]);
    G = createpoints(G,PD,clD,3:6);
    if ischarin('r',varargin)
        % rectangular notch
        P{1} = PC{2} + [0,c/2];
        P{2} = PC{2} - [0,c/2];
        G = createpoints(G,P,clC,[8,1]);
        G = createcontour(G,1:8,1:8,1);
        G = createphysicalline(G,[1 7 8],1);
    elseif ischarin('v',varargin)
        % V (triangular) notch
        P = PC{2};
        G = createpoints(G,P,clC,1);
        G = createcontour(G,1:7,1:7,1);
        G = createphysicalline(G,[1 7],1);
    else%if ischarin('c',varargin)
        % circular (rounded) notch
        P{1} = PC{2} + [-c/2,c/2];
        P{2} = PC{2};
        P{3} = PC{2} - [c/2,c/2];
        P{4} = PC{2} - [c/2,0];
        G = createpoints(G,P,clC,[8,9,1,10]);
        G = createcircle(G,10,8:9,1);
        G = createcircle(G,10,[9,1],2);
        G = createlines(G,[1:7;2:8]',3:9);
        G = createlineloop(G,1:9,1);
        G = createphysicalline(G,[1 2 3 9],1);
    end
    G = createplanesurface(G,1,1);
    if ischarin('recombine',varargin)
        G = recombinesurface(G,1);
    end
    G = createphysicalsurface(G,1,1);
    
elseif indim==3
    PB{1} = PC{1} + [0,c/2,0];
    PB{2} = PC{1} - [0,c/2,0];
    PB{3} = PC{4} + [0,c/2,0];
    PB{4} = PC{4} - [0,c/2,0];
    G = GMSHFILE();
    if ischarin('refinecrack',varargin)
        clB = clC;
    else
        clB = clD;
    end
    G = createpoints(G,PB,clB,[6 4 5 3]);
    if ischarin('r',varargin)
        % rectangular cuboid notch
        P{1} = PC{2} + [0,c/2,0];
        P{2} = PC{2} - [0,c/2,0];
        P{3} = PC{3} + [0,c/2,0];
        P{4} = PC{3} - [0,c/2,0];
        
        G = createpoints(G,PD,clD,9:16);
        G = createpoints(G,P,clC,[7 1 8 2]);
        G = createcontour(G,[2 1 4 3],1:4,1);
        G = createplanesurface(G,1,1);
        
        G = createcontour(G,5:8,5:8,2);
        G = createplanesurface(G,2,2);
        
        G = createlines(G,[[7 1];[2 8]],9:10);
        G = createlineloop(G,[9 -1 10 -7],3);
        G = createplanesurface(G,3,3);
        
        G = createlines(G,[[3 13];[13 14];[14 15];[15 16];[16 5]],11:15);
        G = createlineloop(G,[-4 11:15 -8 -10],4);
        G = createplanesurface(G,4,4);
        
        G = createlines(G,[[6 12];[12 11];[11 10];[10 9];[9 4]],16:20);
        G = createlineloop(G,[-6 16:20 -2 -9],5);
        G = createplanesurface(G,5,5);
        
        G = createlines(G,[[14 10];[9 13]],21:22);
        G = createlineloop(G,-[12 21 19 22],6);
        G = createplanesurface(G,6,6);
        
        G = createlines(G,[[15 11];[12 16]],23:24);
        G = createlineloop(G,[-14 23 -17 24],7);
        G = createplanesurface(G,7,7);
        
        G = createlineloop(G,[-15 -24 -16 -5],8);
        G = createplanesurface(G,8,8);
        
        G = createlineloop(G,[-11 -3 -20 22],9);
        G = createplanesurface(G,9,9);
        
        G = createlineloop(G,[21 -18 -23 -13],10);
        G = createplanesurface(G,10,10);
        
        if ischarin('recombine',varargin)
            G = recombinesurface(G,1);
            G = recombinesurface(G,2);
            G = recombinesurface(G,3);
            G = recombinesurface(G,4);
            G = recombinesurface(G,5);
            G = recombinesurface(G,6);
            G = recombinesurface(G,7);
            G = recombinesurface(G,8);
            G = recombinesurface(G,9);
            G = recombinesurface(G,10);
        end
        G = createsurfaceloop(G,1:10,1);
        G = createphysicalsurface(G,[1 2 3],1);

    elseif ischarin('v',varargin)
        % V (triangular) notch
        G = createpoints(G,PD,clD,7:14);
        G = createpoints(G,PC(2:3),clC,1:2);
        G = createcontour(G,[2 1 4 3],1:4,1);
        G = createplanesurface(G,1,1);
        
        G = createlines(G,[[2 5];[5 6];[6 1]],5:7);
        G = createlineloop(G,[5:7 -1],2);
        G = createplanesurface(G,2,2);
        
        G = createlines(G,[[3 11];[11 12];[12 13];[13 14];[14 5]],8:12);
        G = createlineloop(G,[-4 8:12 -5],3);
        G = createplanesurface(G,3,3);
        
        G = createlines(G,[[6 10];[10 9];[9 8];[8 7];[7 4]],13:17);
        G = createlineloop(G,[-7 13:17 -2],4);
        G = createplanesurface(G,4,4);
        
        G = createlines(G,[[12 8];[7 11]],18:19);
        G = createlineloop(G,-[9 18 16 19],5);
        G = createplanesurface(G,5,5);
        
        G = createlines(G,[[13 9];[10 14]],20:21);
        G = createlineloop(G,[-11 20 -14 21],6);
        G = createplanesurface(G,6,6);
        
        G = createlineloop(G,[-12 -21 -13 -6],7);
        G = createplanesurface(G,7,7);
        
        G = createlineloop(G,[-8 -3 -17 19],8);
        G = createplanesurface(G,8,8);
        
        G = createlineloop(G,[18 -15 -20 -10],9);
        G = createplanesurface(G,9,9);
        
        if ischarin('recombine',varargin)
            G = recombinesurface(G,1);
            G = recombinesurface(G,2);
            G = recombinesurface(G,3);
            G = recombinesurface(G,4);
            G = recombinesurface(G,5);
            G = recombinesurface(G,6);
            G = recombinesurface(G,7);
            G = recombinesurface(G,8);
            G = recombinesurface(G,9);
        end
        G = createsurfaceloop(G,1:9,1);
        G = createphysicalline(G,1,1);
        G = createphysicalsurface(G,[1 2],1);

    else%if ischarin('c',varargin)
        % circular (rounded) cuboid notch
        P{1} = PC{2} + [-c/2,c/2,0];
        P{2} = PC{2};
        P{3} = PC{2} - [c/2,c/2,0];
        P{4} = PC{2} - [c/2,0,0];
        P{5} = PC{3} + [-c/2,c/2,0];
        P{6} = PC{3};
        P{7} = PC{3} - [c/2,c/2,0];
        P{8} = PC{3} - [c/2,0,0];
        
        G = createpoints(G,PD,clD,9:16);
        G = createpoints(G,P,clC,[7 17 1 18 8 19 2 20]);
        G = createcontour(G,[2 1 4 3],1:4,1);
        G = createplanesurface(G,1,1);
        
        G = createcontour(G,5:8,5:8,2);
        G = createplanesurface(G,2,2);
        
        G = createline(G,[17 19],10);
        G = createcircle(G,18,[7 17],9);
        G = createcircle(G,20,[8 19],25);
        G = createcurveloop(G,[9 10 -25 -7],3);
        G = createsurface(G,3,3);
        
        G = createcircle(G,18,[17 1],26);
        G = createcircle(G,20,[19 2],27);
        G = createcurveloop(G,[26 -1 -27 -10],4);
        G = createsurface(G,4,4);
        
        G = createlines(G,[[3 13];[13 14];[14 15];[15 16];[16 5]],11:15);
        G = createlineloop(G,[-4 11:15 -8 25 27],5);
        G = createplanesurface(G,5,5);
        
        G = createlines(G,[[6 12];[12 11];[11 10];[10 9];[9 4]],16:20);
        G = createlineloop(G,[-6 16:20 -2 -26 -9],6);
        G = createplanesurface(G,6,6);
        
        G = createlines(G,[[14 10];[9 13]],21:22);
        G = createlineloop(G,-[12 21 19 22],7);
        G = createplanesurface(G,7,7);
        
        G = createlines(G,[[15 11];[12 16]],23:24);
        G = createlineloop(G,[-14 23 -17 24],8);
        G = createplanesurface(G,8,8);
        
        G = createlineloop(G,[-15 -24 -16 -5],9);
        G = createplanesurface(G,9,9);
        
        G = createlineloop(G,[-11 -3 -20 22],10);
        G = createplanesurface(G,10,10);
        
        G = createlineloop(G,[21 -18 -23 -13],11);
        G = createplanesurface(G,11,11);
        
        if ischarin('recombine',varargin)
            G = recombinesurface(G,1);
            G = recombinesurface(G,2);
            G = recombinesurface(G,3);
            G = recombinesurface(G,4);
            G = recombinesurface(G,5);
            G = recombinesurface(G,6);
            G = recombinesurface(G,7);
            G = recombinesurface(G,8);
            G = recombinesurface(G,9);
            G = recombinesurface(G,10);
            G = recombinesurface(G,11);
        end
        G = createsurfaceloop(G,1:11,1);
        G = createphysicalsurface(G,[1 2 3 4],1);

    end
    G = createvolume(G,1,1);
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
