function testportiquejc(r,varargin)

L = 1;
H = 1;
P = POINT([0,0,0;L,0,0;L,L,0;0,L,0;0,0,H;L,0,H;L,L,H;0,L,H]);

clear LI;
LI{1} = mesh(LIGNE(P(1),P(5)),r);
LI{2} = mesh(LIGNE(P(2),P(6)),r);
LI{3} = mesh(LIGNE(P(3),P(7)),r);
LI{4} = mesh(LIGNE(P(4),P(8)),r);
LI{5} = mesh(LIGNE(P(5),P(6)),r);
LI{6} = mesh(LIGNE(P(6),P(7)),r);
LI{7} = mesh(LIGNE(P(7),P(8)),r);
LI{8} = mesh(LIGNE(P(8),P(5)),r);

LI = union(LI{:});
ray = 0.1;
mat = ELAS_BEAM('RHO',7600,'E',70e9,'NU',0,'S',pi*ray^2,'IY',pi*ray^4/4,'IZ',pi*ray^4/4,'IX',pi*ray^4/2);

S = convertelem(LI,'BEAM',mat);
S = final(S);

PL1 = PLAN(P(1),P(2),P(3));
S = addcl(S,PL1,{'U','R'},0);

K = calc_rigi(S);
M = calc_mass(S);

m = getcharin('mode',varargin,[]);
if m
    [V,D]=calc_mode(K,M,m);
    axis0 =[-L/3,L+L/3,-H/3,H+H/3];
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
    
    f = nodalload(S,P([5:8]'),'FY',1e4);

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
