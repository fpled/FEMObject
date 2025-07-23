function varargout=testpout3D(r,varargin)
close all
L=1;
d=0.1;

E=70e9;
I=pi*d^4/64;
Sec = pi*d^2/4;
mat = ELAS_BEAM('E',E,'NU',0.3,'RHO',7600,'S',Sec,'IZ',I,'IY',I,'IX',I*2);

P1=POINT([0,0,0]);
P2=POINT([L,0,0]);
P3=POINT([L/2,0,0]);
S1=mesh(LINE(P1,P2),r,'material',mat);
%S=MODEL('TRID');
%S=addelem(S,'BEAM',S1,mat,'param',VECTOR([0;1;0]));
S=convertelem(S1,'BEAM','param',VECTOR([0;1;0]));
S=final(S,'norenum')
S=addcl(S,P1,{'U','R'},0);
%S=addcl(S,P2,{'U','R'},1);

K=calc_rigi(S);

m = getcharin('mode',varargin,[]);
if m
    ampl=getcharin('ampl',varargin,0.1);
    M=calc_mass(S);
%     m=varargin{pos+1};
    [V,D]=calc_mode(K,M,m);
    n=10; ampli = cos(linspace(0,2*pi,n))*ampl;
    figure(1);
    axis0 =[-L/5,L+L/5,-L/2,L/2];
    varargout{1}=V;
    varargout{2}=D;

    for i=1:length(ampli) 
        clf;
        for k=m(:)'
            subplot(1,length(m),find(k==m))
            plot(S+ampli(i)*V(:,k));

            axis(axis0);
        end
%FILM(i)=getframe;
        pause(0.01);
    end
else
% movie(FILM,1,24)

    ampl=getcharin('ampl',varargin,1e5);

    F = 1;
    f=nodalload(S,P2,'FY',F);    
    q=K\f;
    q=unfreevector(S,q);
    figure(2)
    clf
    varargout{1}=q;
    plot(S,'EdgeColor','b')
    plot(S+ampl*q)

    a=F*L^3/3/getparam(mat,'E')/getparam(mat,'IZ');
    b=q(findddl(S,'UY',P2))

    s=calc_sigma(S,q,'smooth');
    fprintf('systeme de coordonnees de la section\n')
    getsyscoordlocal(S.groupelem{1},1)

    figure(6)
    clf
    subplot(1,4,1)
    plot(s,S,'compo','EFFX')
    colorbar('horiz')
    axis image
    subplot(1,4,2)
    plot(s,S,'compo','MOMX')
    colorbar('horiz')
    subplot(1,4,3)
    plot(s,S,'compo','MOMY')
    colorbar('horiz')
    subplot(1,4,4)
    plot(s,S,'compo','MOMZ')
    colorbar('horiz')

    fprintf('error = %d\n',abs((a-b)/b))

end
%s=calc_sigma(q,M);
