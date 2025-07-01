function varargout = gmshDomainWithTwoAsymmetricEdgeNotches(D,Ca,Cb,c,clD,clC,filename,indim,varargin)
% function varargout = gmshDomainWithTwoAsymmetricEdgeNotches(D,Ca,Cb,c,clD,clC,filename,indim,notchtype)
% D : DOMAIN
% Ca, Cb : LIGNE in dim 2, QUADRANGLE in dim 3
% c : width of the edge notch
% clD, clC : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getindim(D) by default)
% notchtype : 'c', 'circ' or 'circular' for circular notch (by default)
%             'r', 'rect' or 'rectangular' for rectangular notch
%             'v', 'V' or 'triangular' for V notch

if nargin<8 || isempty(indim)
    indim = getindim(D);
end
if nargin<6 || isempty(clC)
    clC = clD;
end

PD = getvertices(D);
PCa = getvertices(Ca);
PCb = getvertices(Cb);

if indim==2
    PCa1 = min(PCa{:});
    PCb2 = max(PCb{:});
    PB{1} = PCa1 + [0,c/2];
    PB{2} = PCa1 - [0,c/2];
    PB{3} = PCb2 - [0,c/2];
    PB{4} = PCb2 + [0,c/2];
    G = GMSHFILE();
    if ischarin('refinecrack',varargin)
        clB = clC;
    else
        clB = clD;
    end
    PCa2 = max(PCa{:});
    PCb1 = min(PCb{:});
    if ischarin('r',varargin) || ischarin('rect',varargin) || ischarin('rectangular',varargin)
        % rectangular notch
        G = createpoints(G,PB,clB,[11,2,5,8]);
        G = createpoints(G,PD,clD,[3:4,9:10]);
        P{1} = PCa2 + [0,c/2];
        P{2} = PCa2 - [0,c/2];
        P{3} = PCb1 - [0,c/2];
        P{4} = PCb1 + [0,c/2];
        G = createpoints(G,P,clC,[12,1,6,7]);
        G = createcontour(G,1:12,1:12,1);
        % G = createphysicalcurve(G,[1 5:7 11 12],1);
    elseif ischarin('v',varargin) || ischarin('V',varargin)|| ischarin('triangular',varargin)
        % V (triangular) notch
        G = createpoints(G,PB,clB,[10,2,5,7]);
        G = createpoints(G,PD,clD,[3:4,8:9]);
        P{1} = PCa2;
        P{2} = PCb1;
        G = createpoints(G,P,clC,[1,6]);
        G = createcontour(G,1:10,1:10,1);
        % G = createphysicalcurve(G,[1 5 6 10],1);
    else%if ischarin('c',varargin) || ischarin('circ',varargin) || ischarin('circular',varargin)
        % circular (rounded) notch
        G = createpoints(G,PB,clB,[12,2,5,9]);
        G = createpoints(G,PD,clD,[3:4,10:11]);
        P{1} = PCa2 + [-c/2,c/2];
        P{2} = PCa2;
        P{3} = PCa2 - [c/2,c/2];
        P{4} = PCa2 - [c/2,0];
        P{5} = PCb1 + [c/2,-c/2];
        P{6} = PCb1;
        P{7} = PCb1 + [c/2,c/2];
        P{8} = PCb1 + [c/2,0];
        G = createpoints(G,P,clC,[13,14,1,15,6,7,8,16]);
        G = createcircle(G,15,13:14,1);
        G = createcircle(G,15,[14,1],2);
        G = createcircle(G,16,6:7,8);
        G = createcircle(G,16,7:8,9);
        G = createlines(G,[1:5 8:12;2:6 9:13]',[3:7 10:14]);
        G = createcurveloop(G,1:14,1);
        G = createphysicalcurve(G,[1:3 7:10 14],1);
    end
    G = createplanesurface(G,1,1);
    if ischarin('recombine',varargin)
        G = recombinesurface(G,1);
    end
    G = createphysicalsurface(G,1,1);
    
elseif indim==3
    PB{1} = PCa{1} + [0,c/2,0];
    PB{2} = PCa{1} - [0,c/2,0];
    PB{3} = PCa{4} + [0,c/2,0];
    PB{4} = PCa{4} - [0,c/2,0];
    PB{5} = PCb{1} + [0,c/2,0];
    PB{6} = PCb{1} - [0,c/2,0];
    PB{7} = PCb{4} + [0,c/2,0];
    PB{8} = PCb{4} - [0,c/2,0];
    G = GMSHFILE();
    if ischarin('refinecrack',varargin)
        clB = clC;
    else
        clB = clD;
    end
    if ischarin('r',varargin) || ischarin('rect',varargin) || ischarin('rectangular',varargin)
        % rectangular cuboid notch
        G = createpoints(G,PB,clB,[6 4 5 3 14 12 13 11]);
        P{1} = PCa{2} + [0,c/2,0];
        P{2} = PCa{2} - [0,c/2,0];
        P{3} = PCa{3} + [0,c/2,0];
        P{4} = PCa{3} - [0,c/2,0];
        P{5} = PCb{2} + [0,c/2,0];
        P{6} = PCb{2} - [0,c/2,0];
        P{7} = PCb{3} + [0,c/2,0];
        P{8} = PCb{3} - [0,c/2,0];
        
        G = createpoints(G,PD,clD,17:24);
        G = createpoints(G,P,clC,[7 1 8 2 15 9 16 10]);
        G = createcontour(G,[2 1 4 3],1:4,1);
        G = createplanesurface(G,1,1);
        
        G = createcontour(G,5:8,5:8,2);
        G = createplanesurface(G,2,2);
        
        G = createlines(G,[[7 1];[2 8]],9:10);
        G = createcurveloop(G,[9 -1 10 -7],3);
        G = createplanesurface(G,3,3);
        
        G = createcontour(G,9:12,11:14,4);
        G = createplanesurface(G,4,4);
        
        G = createcontour(G,[14 13 16 15],15:18,5);
        G = createplanesurface(G,5,5);
        
        G = createlines(G,[[15 9];[10 16]],19:20);
        G = createcurveloop(G,-[19 17 20 11],6);
        G = createplanesurface(G,6,6);

        G = createlines(G,[[3 21];[21 22];[22 11];[13 23];[23 24];[24 5]],21:26);
        G = createcurveloop(G,[-4 21:23 -12 20 -16 24:26 -8 -10],7);
        G = createplanesurface(G,7,7);
        
        G = createlines(G,[[6 20];[20 19];[19 14];[12 18];[18 17];[17 4]],27:32);
        G = createcurveloop(G,[-6 27:29 -18 19 -14 30:32 -2 -9],8);
        G = createplanesurface(G,8,8);
        
        G = createlines(G,[[22 18];[17 21]],33:34);
        G = createcurveloop(G,-[22 34 31 33],9);
        G = createplanesurface(G,9,9);
        
        G = createlines(G,[[23 19];[20 24]],35:36);
        G = createcurveloop(G,[-25 35 -28 36],10);
        G = createplanesurface(G,10,10);
        
        G = createcurveloop(G,[-26 -36 -27 -5],11);
        G = createplanesurface(G,11,11);
        
        G = createcurveloop(G,[-21 -3 -32 34],12);
        G = createplanesurface(G,12,12);
        
        G = createcurveloop(G,[-15 -29 -35 -24],13);
        G = createplanesurface(G,13,13);
        
        G = createcurveloop(G,[33 -30 -13 -23],14);
        G = createplanesurface(G,14,14);
        
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
            G = recombinesurface(G,12);
            G = recombinesurface(G,13);
            G = recombinesurface(G,14);
        end
        G = createsurfaceloop(G,1:14,1);
        % G = createphysicalsurface(G,[1 2 3 4 5 6],1);

    elseif ischarin('v',varargin) || ischarin('V',varargin)|| ischarin('triangular',varargin)
        % V (triangular) notch
        G = createpoints(G,PB,clB,[6 4 5 3 12 10 11 9]);
        G = createpoints(G,PD,clD,13:20);
        G = createpoints(G,[PCa(2:3),PCb(2:3)],clC,[1 2 7 8]);
        G = createcontour(G,[2 1 4 3],1:4,1);
        G = createplanesurface(G,1,1);
        
        G = createlines(G,[[2 5];[5 6];[6 1]],5:7);
        G = createcurveloop(G,[5:7 -1],2);
        G = createplanesurface(G,2,2);
        
        G = createcontour(G,7:10,8:11,3);
        G = createplanesurface(G,3,3);
        
        G = createlines(G,[[7 12];[12 11];[11 8]],12:14);
        G = createcurveloop(G,[12:14 -8],4);
        G = createplanesurface(G,4,4);
        
        G = createlines(G,[[3 17];[17 18];[18 9];[11 19];[19 20];[20 5]],15:20);
        G = createcurveloop(G,[-4 15:17 -9 -14 18:20 -5],5);
        G = createplanesurface(G,5,5);
        
        G = createlines(G,[[6 16];[16 15];[15 12];[10 14];[14 13];[13 4]],21:26);
        G = createcurveloop(G,[-7 21:23 -12 -11 24:26 -2],6);
        G = createplanesurface(G,6,6);
        
        G = createlines(G,[[18 14];[13 17]],27:28);
        G = createcurveloop(G,-[16 28 25 27],7);
        G = createplanesurface(G,7,7);
        
        G = createlines(G,[[19 15];[16 20]],29:30);
        G = createcurveloop(G,[-19 29 -22 30],8);
        G = createplanesurface(G,8,8);
        
        G = createcurveloop(G,[-20 -30 -21 -6],9);
        G = createplanesurface(G,9,9);
        
        G = createcurveloop(G,[-15 -3 -26 28],10);
        G = createplanesurface(G,10,10);
        
        G = createcurveloop(G,[-13 -23 -29 -18],11);
        G = createplanesurface(G,11,11);
        
        G = createcurveloop(G,[27 -24 -10 -17],12);
        G = createplanesurface(G,12,12);
        
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
            G = recombinesurface(G,12);
        end
        G = createsurfaceloop(G,1:12,1);
        % G = createphysicalcurve(G,[1 8],1);
        % G = createphysicalsurface(G,[1 2 3 4],1);

    else%if ischarin('c',varargin) || ischarin('circ',varargin) || ischarin('circular',varargin)
        % circular (rounded) cuboid notch
        G = createpoints(G,PB,clB,[6 4 5 3 14 12 13 11]);
        P{1} = PCa{2} + [-c/2,c/2,0];
        P{2} = PCa{2};
        P{3} = PCa{2} - [c/2,c/2,0];
        P{4} = PCa{2} - [c/2,0,0];
        P{5} = PCa{3} + [-c/2,c/2,0];
        P{6} = PCa{3};
        P{7} = PCa{3} - [c/2,c/2,0];
        P{8} = PCa{3} - [c/2,0,0];
        P{9} = PCb{2} + [c/2,c/2,0];
        P{10} = PCb{2};
        P{11} = PCb{2} + [c/2,-c/2,0];
        P{12} = PCb{2} + [c/2,0,0];
        P{13} = PCb{3} + [c/2,c/2,0];
        P{14} = PCb{3};
        P{15} = PCb{3} + [c/2,-c/2,0];
        P{16} = PCb{3} + [c/2,0,0];
        
        G = createpoints(G,PD,clD,17:24);
        G = createpoints(G,P,clC,[7 25 1 26 8 27 2 28 15 29 9 30 16 31 10 32]);
        G = createcontour(G,[2 1 4 3],1:4,1);
        G = createplanesurface(G,1,1);
        
        G = createcontour(G,5:8,5:8,2);
        G = createplanesurface(G,2,2);
        
        G = createline(G,[25 27],10);
        G = createcircle(G,26,[7 25],9);
        G = createcircle(G,28,[8 27],37);
        G = createcurveloop(G,[9 10 -37 -7],3);
        G = createsurface(G,3,3);
        
        G = createcircle(G,26,[25 1],38);
        G = createcircle(G,28,[27 2],39);
        G = createcurveloop(G,[38 -1 -39 -10],4);
        G = createsurface(G,4,4);
        
        G = createcontour(G,9:12,11:14,5);
        G = createplanesurface(G,5,5);
        
        G = createcontour(G,[14 13 16 15],15:18,6);
        G = createplanesurface(G,6,6);
        
        G = createline(G,[29 31],20);
        G = createcircle(G,30,[15 29],19);
        G = createcircle(G,32,[16 31],40);
        G = createcurveloop(G,[-19 -17 40 -20],7);
        G = createsurface(G,7,7);
        
        G = createcircle(G,30,[29 9],41);
        G = createcircle(G,32,[31 10],42);
        G = createcurveloop(G,[-41 20 42 -11],8);
        G = createsurface(G,8,8);
        
        G = createlines(G,[[3 21];[21 22];[22 11];[13 23];[23 24];[24 5]],21:26);
        G = createcurveloop(G,[-4 21:23 -12 -42 -40 -16 24:26 -8 37 39],9);
        G = createplanesurface(G,9,9);
        
        G = createlines(G,[[6 20];[20 19];[19 14];[12 18];[18 17];[17 4]],27:32);
        G = createcurveloop(G,[-6 27:29 -18 19 41 -14 30:32 -2 -38 -9],10);
        G = createplanesurface(G,10,10);
        
        G = createlines(G,[[22 18];[17 21]],33:34);
        G = createcurveloop(G,-[22 33 31 34],11);
        G = createplanesurface(G,11,11);
        
        G = createlines(G,[[23 19];[20 24]],35:36);
        G = createcurveloop(G,[-25 35 -28 36],12);
        G = createplanesurface(G,12,12);
        
        G = createcurveloop(G,[-26 -36 -27 -5],13);
        G = createplanesurface(G,13,13);
        
        G = createcurveloop(G,[-21 -3 -32 34],14);
        G = createplanesurface(G,14,14);
        
        G = createcurveloop(G,[-15 -29 -35 -24],15);
        G = createplanesurface(G,15,15);
        
        G = createcurveloop(G,[33 -30 -13 -23],16);
        G = createplanesurface(G,16,16);
        
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
            G = recombinesurface(G,12);
            G = recombinesurface(G,13);
            G = recombinesurface(G,14);
            G = recombinesurface(G,15);
            G = recombinesurface(G,16);
        end
        G = createsurfaceloop(G,1:16,1);
        % G = createphysicalsurface(G,[1 2 3 4 5 6 7 8],1);

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

if nargin>=7 && ischar(filename)
    G = setfile(G,filename);
end

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(D):-1:getdim(D)-n+1,varargin{:});
