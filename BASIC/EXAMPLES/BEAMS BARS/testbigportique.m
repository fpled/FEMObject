function testbigportique(r,varargin)

L = 1;
H = 1;
nblock = 5;
HPC = 3*(H*nblock);
PC = POINT([L/2,L/2,HPC]);
a = 3;
% V = VECTEUR([1,1,-1,-1;1,-1,-1,1;-a,-a,-a,-a]);
Pfloor = POINT([0,0,0;L,0,0;L,L,0;0,L,0]);
V = (PC-Pfloor)*(H/HPC);
for k=0:nblock
    Pstage{k+1} = Pfloor+V*k;    
end
Ptop = Pstage{end};
clear LI;
for k=1:nblock
    P1 = Pstage{k};
    P2 = Pstage{k+1};
    P = [P1,P2];
    LI{1} = mesh(LIGNE(P(1),P(5)),r);
    LI{2} = mesh(LIGNE(P(2),P(6)),r);
    LI{3} = mesh(LIGNE(P(3),P(7)),r);
    LI{4} = mesh(LIGNE(P(4),P(8)),r);
    LI{5} = mesh(LIGNE(P(5),P(6)),r);
    LI{6} = mesh(LIGNE(P(6),P(7)),r);
    LI{7} = mesh(LIGNE(P(7),P(8)),r);
    LI{8} = mesh(LIGNE(P(8),P(5)),r);
    LI{9} = mesh(LIGNE(P(1),P(7)),r);
    LI{10} = mesh(LIGNE(P(2),P(8)),r);
    LI{11} = mesh(LIGNE(P(3),P(5)),r);
    LI{12} = mesh(LIGNE(P(4),P(6)),r);
    if k==1
        S = union(LI{:});
    else
        S = union(S,union(LI{:})); 
    end
end

ray = 0.1;
mat = MATERIALS();
mat{1} = ELAS_BEAM('RHO',7600,'E',70e9,'NU',0,'S',pi*ray^2,'IY',pi*ray^4/4,'IZ',pi*ray^4/4,'IX',pi*ray^4/2);

S = concatgroupelem(S);
S = convertelem(S,'BEAM',mat{1});
S = final(S);

PL1  = PLAN(Pfloor(1),Pfloor(2),Pfloor(3));
S = addcl(S,PL1,{'U','R'},0);

K = calc_rigi(S);
M = calc_mass(S);

m = getcharin('mode',varargin,[]);
if m
    [V,D] = calc_mode(K,M,m);
    axis0 = [-L/3,L+L/3,-H/3,H+H/3];
    ampl = getcharin('ampl',varargin,5);
    n = 30;
    ampli = cos(linspace(0,2*pi,n))*ampl;
    figure(1)
    for i=1:length(ampli) 
        clf
        for j=1:length(m)
            subplot(ceil(sqrt(length(m))),ceil(sqrt(length(m))),j);
            plot(S+ampli(i)*V(:,j));
            axis(axis0);
        end
        pause(1/n);
    end
else
    ampl = getcharin('ampl',varargin,1e4);

    f = nodalload(S,Ptop,'FY',1e4);

    q = K\f;
    
    figure(2)
    clf
    plot(S,'EdgeColor','b','node')
    plot(S+ampl*q,'node');

    s = calc_sigma(S,q,'node');

    figure(3)
    clf
    subplot(1,3,1)
    plot(s,S+ampl*q,'compo','MOMX','node')
    colorbar('h')
    subplot(1,3,2)
    plot(s,S+ampl*q,'compo','MOMY','node')
    colorbar('h')
    subplot(1,3,3)
    plot(s,S+ampl*q,'compo','MOMZ','node')
    colorbar('h')
    
    figure(4)
    clf
    plot(s,S+ampl*q,'compo','EFFX','node')
    colorbar
end
