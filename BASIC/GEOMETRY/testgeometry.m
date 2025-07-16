center = [0,0,0];

vx = [1,0];
vy = [0,1];
vxy = [1,1];

nx = [1,0,0];
ny = [0,1,0];
nz = [0,0,1];
nzm = [0,0,-1];
nxyz = [1,1,1];

%% CIRCLE
r = 0.3;
ampl = r;

C = CIRCLE(center(1),center(2),center(3),r);

Cvy = CIRCLE(center(1),center(2),center(3),r,vy(1),vy(2));
Cvxy = CIRCLE(center(1),center(2),center(3),r,vxy(1),vxy(2));

Cnx = CIRCLE(center(1),center(2),center(3),r,nx(1),nx(2),nx(3));
Cnxvy = CIRCLE(center(1),center(2),center(3),r,nx(1),nx(2),nx(3),vy(1),vy(2));
Cnxvxy = CIRCLE(center(1),center(2),center(3),r,nx(1),nx(2),nx(3),vxy(1),vxy(2));

Cny = CIRCLE(center(1),center(2),center(3),r,ny(1),ny(2),ny(3));
Cnyvy = CIRCLE(center(1),center(2),center(3),r,ny(1),ny(2),ny(3),vy(1),vy(2));
Cnyvxy = CIRCLE(center(1),center(2),center(3),r,ny(1),ny(2),ny(3),vxy(1),vxy(2));

Cnxyz = CIRCLE(center(1),center(2),center(3),r,nxyz(1),nxyz(2),nxyz(3));
Cnxyzvy = CIRCLE(center(1),center(2),center(3),r,nxyz(1),nxyz(2),nxyz(3),vy(1),vy(2));
Cnxyzvxy = CIRCLE(center(1),center(2),center(3),r,nxyz(1),nxyz(2),nxyz(3),vxy(1),vxy(2));

Cnzm = CIRCLE(center(1),center(2),center(3),r,nzm(1),nzm(2),nzm(3));
Cnzmvy = CIRCLE(center(1),center(2),center(3),r,nzm(1),nzm(2),nzm(3),vy(1),vy(2));
Cnzmvxy = CIRCLE(center(1),center(2),center(3),r,nzm(1),nzm(2),nzm(3),vxy(1),vxy(2));

figure('Name','Circle')
plot(C,'EdgeColor','k')
hold on
plot(getnormal(C),'k',ampl)
% plot(Cvy,'EdgeColor','k')
% plot(Cvxy,'EdgeColor','k')

plot(Cnx,'EdgeColor','b')
plot(getnormal(Cnx),'b',ampl)
% plot(Cnxvy,'EdgeColor','b')
% plot(Cnxvxy,'EdgeColor','b')

plot(Cny,'EdgeColor','r')
plot(getnormal(Cny),'r',ampl)
% plot(Cnyvy,'EdgeColor','r')
% plot(Cnyvxy,'EdgeColor','r')

plot(Cnxyz,'EdgeColor','g')
plot(getnormal(Cnxyz),'g',ampl)
% plot(Cnxyzvy,'EdgeColor','g')
% plot(Cnxyzvxy,'EdgeColor','g')

plot(Cnzm,'EdgeColor','c')
plot(getnormal(Cnzm),'c',ampl)
% plot(Cnzmvy,'EdgeColor','c')
% plot(Cnzmvxy,'EdgeColor','c')

hold off
grid on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')

%% ELLIPSE
a = 0.5;
b = 0.2;
ampl = max(a,b);

E = ELLIPSE(center(1),center(2),center(3),a,b);

Evy = ELLIPSE(center(1),center(2),center(3),a,b,vy(1),vy(2));
Evxy = ELLIPSE(center(1),center(2),center(3),a,b,vxy(1),vxy(2));

Enx = ELLIPSE(center(1),center(2),center(3),a,b,nx(1),nx(2),nx(3));
Enxvy = ELLIPSE(center(1),center(2),center(3),a,b,nx(1),nx(2),nx(3),vy(1),vy(2));
Enxvxy = ELLIPSE(center(1),center(2),center(3),a,b,nx(1),nx(2),nx(3),vxy(1),vxy(2));

Eny = ELLIPSE(center(1),center(2),center(3),a,b,ny(1),ny(2),ny(3));
Enyvy = ELLIPSE(center(1),center(2),center(3),a,b,ny(1),ny(2),ny(3),vy(1),vy(2));
Enyvxy = ELLIPSE(center(1),center(2),center(3),a,b,ny(1),ny(2),ny(3),vxy(1),vxy(2));

Enxyz = ELLIPSE(center(1),center(2),center(3),a,b,nxyz(1),nxyz(2),nxyz(3));
Enxyzvy = ELLIPSE(center(1),center(2),center(3),a,b,nxyz(1),nxyz(2),nxyz(3),vy(1),vy(2));
Enxyzvxy = ELLIPSE(center(1),center(2),center(3),a,b,nxyz(1),nxyz(2),nxyz(3),vxy(1),vxy(2));

Enzm = ELLIPSE(center(1),center(2),center(3),a,b,nzm(1),nzm(2),nzm(3));
Enzmvy = ELLIPSE(center(1),center(2),center(3),a,b,nzm(1),nzm(2),nzm(3),vy(1),vy(2));
Enzmvxy = ELLIPSE(center(1),center(2),center(3),a,b,nzm(1),nzm(2),nzm(3),vxy(1),vxy(2));

figure('Name','Ellipse')
plot(E,'EdgeColor','k')
hold on
plot(getnormal(E),'k',ampl)
% plot(Evy,'EdgeColor','k')
% plot(Evxy,'EdgeColor','k')

plot(Enx,'EdgeColor','b')
plot(getnormal(Enx),'b',ampl)
% plot(Enxvy,'EdgeColor','b')
% plot(Enxvxy,'EdgeColor','b')

plot(Eny,'EdgeColor','r')
plot(getnormal(Eny),'r',ampl)
% plot(Enyvy,'EdgeColor','r')
% plot(Enyvxy,'EdgeColor','r')

plot(Enxyz,'EdgeColor','g')
plot(getnormal(Enxyz),'g',ampl)
% plot(Enxyzvy,'EdgeColor','g')
% plot(Enxyzvxy,'EdgeColor','g')

plot(Enzm,'EdgeColor','c')
plot(getnormal(Enzm),'c',ampl)
% plot(Enzmvy,'EdgeColor','c')
% plot(Enzmvxy,'EdgeColor','c')

hold off
grid on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')

%% CYLINDER
r = 0.3;
h = 1;
ampl = h;

CY = CYLINDER(center(1),center(2),center(3),r,h);

CYvy = CYLINDER(center(1),center(2),center(3),r,vy(1),vy(2),h);
CYvxy = CYLINDER(center(1),center(2),center(3),r,vxy(1),vxy(2),h);

CYnx = CYLINDER(center(1),center(2),center(3),r,nx(1),nx(2),nx(3),h);
CYnxvy = CYLINDER(center(1),center(2),center(3),r,nx(1),nx(2),nx(3),vy(1),vy(2),h);
CYnxvxy = CYLINDER(center(1),center(2),center(3),r,nx(1),nx(2),nx(3),vxy(1),vxy(2),h);

CYny = CYLINDER(center(1),center(2),center(3),r,ny(1),ny(2),ny(3),h);
CYnyvy = CYLINDER(center(1),center(2),center(3),r,ny(1),ny(2),ny(3),vy(1),vy(2),h);
CYnyvxy = CYLINDER(center(1),center(2),center(3),r,ny(1),ny(2),ny(3),vxy(1),vxy(2),h);

CYnxyz = CYLINDER(center(1),center(2),center(3),r,nxyz(1),nxyz(2),nxyz(3),h);
CYnxyzvy = CYLINDER(center(1),center(2),center(3),r,nxyz(1),nxyz(2),nxyz(3),vy(1),vy(2),h);
CYnxyzvxy = CYLINDER(center(1),center(2),center(3),r,nxyz(1),nxyz(2),nxyz(3),vxy(1),vxy(2),h);

CYnzm = CYLINDER(center(1),center(2),center(3),r,nzm(1),nzm(2),nzm(3),h);
CYnzmvy = CYLINDER(center(1),center(2),center(3),r,nzm(1),nzm(2),nzm(3),vy(1),vy(2),h);
CYnzmvxy = CYLINDER(center(1),center(2),center(3),r,nzm(1),nzm(2),nzm(3),vxy(1),vxy(2),h);

figure('Name','Cylinder')
plot(CY,'EdgeColor','k')
hold on
plot(getnormal(CY),'k',ampl)
% plot(CYvy,'EdgeColor','k')
% plot(CYvxy,'EdgeColor','k')

plot(CYnx,'EdgeColor','b')
plot(getnormal(CYnx),'b',ampl)
% plot(CYnxvy,'EdgeColor','b')
% plot(CYnxvxy,'EdgeColor','b')

plot(CYny,'EdgeColor','r')
plot(getnormal(CYny),'r',ampl)
% plot(CYnyvy,'EdgeColor','r')
% plot(CYnyvxy,'EdgeColor','r')

plot(CYnxyz,'EdgeColor','g')
plot(getnormal(CYnxyz),'g',ampl)
% plot(CYnxyzvy,'EdgeColor','g')
% plot(CYnxyzvxy,'EdgeColor','g')

plot(CYnzm,'EdgeColor','c')
plot(getnormal(CYnzm),'c',ampl)
% plot(CYnzmvy,'EdgeColor','c')
% plot(CYnzmvxy,'EdgeColor','c')

hold off
grid on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')

figure('Name','Cylinder')
surf(CY,'EdgeColor','k')
hold on
plot(getnormal(CY),'k',ampl)
% surf(CYvy,'EdgeColor','k')
% surf(CYvxy,'EdgeColor','k')

surf(CYnx,'EdgeColor','b')
plot(getnormal(CYnx),'b',ampl)
% surf(CYnxvy,'EdgeColor','b')
% surf(CYnxvxy,'EdgeColor','b')

surf(CYny,'EdgeColor','r')
plot(getnormal(CYny),'r',ampl)
% surf(CYnyvy,'EdgeColor','r')
% surf(CYnyvxy,'EdgeColor','r')

surf(CYnxyz,'EdgeColor','g')
plot(getnormal(CYnxyz),'g',ampl)
% surf(CYnxyzvy,'EdgeColor','g')
% surf(CYnxyzvxy,'EdgeColor','g')

surf(CYnzm,'EdgeColor','c')
plot(getnormal(CYnzm),'c',ampl)
% surf(CYnzmvy,'EdgeColor','c')
% surf(CYnzmvxy,'EdgeColor','c')

hold off
grid on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')

%% SPHERE
r = 0.3;
ampl = r;

S = SPHERE(center(1),center(2),center(3),r);

Svy = SPHERE(center(1),center(2),center(3),r,vy(1),vy(2));
Svxy = SPHERE(center(1),center(2),center(3),r,vxy(1),vxy(2));

Snx = SPHERE(center(1),center(2),center(3),r,nx(1),nx(2),nx(3));
Snxvy = SPHERE(center(1),center(2),center(3),r,nx(1),nx(2),nx(3),vy(1),vy(2));
Snxvxy = SPHERE(center(1),center(2),center(3),r,nx(1),nx(2),nx(3),vxy(1),vxy(2));

Sny = SPHERE(center(1),center(2),center(3),r,ny(1),ny(2),ny(3));
Snyvy = SPHERE(center(1),center(2),center(3),r,ny(1),ny(2),ny(3),vy(1),vy(2));
Snyvxy = SPHERE(center(1),center(2),center(3),r,ny(1),ny(2),ny(3),vxy(1),vxy(2));

Snxyz = SPHERE(center(1),center(2),center(3),r,nxyz(1),nxyz(2),nxyz(3));
Snxyzvy = SPHERE(center(1),center(2),center(3),r,nxyz(1),nxyz(2),nxyz(3),vy(1),vy(2));
Snxyzvxy = SPHERE(center(1),center(2),center(3),r,nxyz(1),nxyz(2),nxyz(3),vxy(1),vxy(2));

Snzm = SPHERE(center(1),center(2),center(3),r,nzm(1),nzm(2),nzm(3));
Snzmvy = SPHERE(center(1),center(2),center(3),r,nzm(1),nzm(2),nzm(3),vy(1),vy(2));
Snzmvxy = SPHERE(center(1),center(2),center(3),r,nzm(1),nzm(2),nzm(3),vxy(1),vxy(2));

figure('Name','Sphere')
plot(S,'EdgeColor','k')
hold on
plot(getnormal(S),'k',ampl)
% plot(Svy,'EdgeColor','k')
% plot(Svxy,'EdgeColor','k')

plot(Snx,'EdgeColor','b')
plot(getnormal(Snx),'b',ampl)
% plot(Snxvy,'EdgeColor','b')
% plot(Snxvxy,'EdgeColor','b')

plot(Sny,'EdgeColor','r')
plot(getnormal(Sny),'r',ampl)
% plot(Snyvy,'EdgeColor','r')
% plot(Snyvxy,'EdgeColor','r')

plot(Snxyz,'EdgeColor','g')
plot(getnormal(Snxyz),'g',ampl)
% plot(Snxyzvy,'EdgeColor','g')
% plot(Snxyzvxy,'EdgeColor','g')

plot(Snzm,'EdgeColor','c')
plot(getnormal(Snzm),'c',ampl)
% plot(Snzmvy,'EdgeColor','c')
% plot(Snzmvxy,'EdgeColor','c')

hold off
grid on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')

figure('Name','Sphere')
surf(S,'EdgeColor','k')
hold on
plot(getnormal(S),'k',ampl)
% surf(Svy,'EdgeColor','k')
% surf(Svxy,'EdgeColor','k')

surf(Snx,'EdgeColor','b')
plot(getnormal(Snx),'b',ampl)
% surf(Snxvy,'EdgeColor','b')
% surf(Snxvxy,'EdgeColor','b')

surf(Sny,'EdgeColor','r')
plot(getnormal(Sny),'r',ampl)
% surf(Snyvy,'EdgeColor','r')
% surf(Snyvxy,'EdgeColor','r')

surf(Snxyz,'EdgeColor','g')
plot(getnormal(Snxyz),'g',ampl)
% surf(Snxyzvy,'EdgeColor','g')
% surf(Snxyzvxy,'EdgeColor','g')

surf(Snzm,'EdgeColor','c')
plot(getnormal(Snzm),'c',ampl)
% surf(Snzmvy,'EdgeColor','c')
% surf(Snzmvxy,'EdgeColor','c')

hold off
grid on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')

%% ELLIPSOID
a = 0.5;
b = 0.2;
c = 0.1;
ampl = c;

E = ELLIPSOID(center(1),center(2),center(3),a,b,c);

Evy = ELLIPSOID(center(1),center(2),center(3),a,b,c,vy(1),vy(2));
Evxy = ELLIPSOID(center(1),center(2),center(3),a,b,c,vxy(1),vxy(2));

Enx = ELLIPSOID(center(1),center(2),center(3),a,b,c,nx(1),nx(2),nx(3));
Enxvy = ELLIPSOID(center(1),center(2),center(3),a,b,c,nx(1),nx(2),nx(3),vy(1),vy(2));
Enxvxy = ELLIPSOID(center(1),center(2),center(3),a,b,c,nx(1),nx(2),nx(3),vxy(1),vxy(2));

Eny = ELLIPSOID(center(1),center(2),center(3),a,b,c,ny(1),ny(2),ny(3));
Enyvy = ELLIPSOID(center(1),center(2),center(3),a,b,c,ny(1),ny(2),ny(3),vy(1),vy(2));
Enyvxy = ELLIPSOID(center(1),center(2),center(3),a,b,c,ny(1),ny(2),ny(3),vxy(1),vxy(2));

Enxyz = ELLIPSOID(center(1),center(2),center(3),a,b,c,nxyz(1),nxyz(2),nxyz(3));
Enxyzvy = ELLIPSOID(center(1),center(2),center(3),a,b,c,nxyz(1),nxyz(2),nxyz(3),vy(1),vy(2));
Enxyzvxy = ELLIPSOID(center(1),center(2),center(3),a,b,c,nxyz(1),nxyz(2),nxyz(3),vxy(1),vxy(2));

Enzm = ELLIPSOID(center(1),center(2),center(3),a,b,c,nzm(1),nzm(2),nzm(3));
Enzmvy = ELLIPSOID(center(1),center(2),center(3),a,b,c,nzm(1),nzm(2),nzm(3),vy(1),vy(2));
Enzmvxy = ELLIPSOID(center(1),center(2),center(3),a,b,c,nzm(1),nzm(2),nzm(3),vxy(1),vxy(2));

figure('Name','Ellipsoid')
plot(E,'EdgeColor','k')
hold on
plot(getnormal(E),'k',ampl)
% plot(Evy,'EdgeColor','k')
% plot(Evxy,'EdgeColor','k')

plot(Enx,'EdgeColor','b')
plot(getnormal(Enx),'b',ampl)
% plot(Enxvy,'EdgeColor','b')
% plot(Enxvxy,'EdgeColor','b')

plot(Eny,'EdgeColor','r')
plot(getnormal(Eny),'r',ampl)
% plot(Enyvy,'EdgeColor','r')
% plot(Enyvxy,'EdgeColor','r')

plot(Enxyz,'EdgeColor','g')
plot(getnormal(Enxyz),'g',ampl)
% plot(Enxyzvy,'EdgeColor','g')
% plot(Enxyzvxy,'EdgeColor','g')

plot(Enzm,'EdgeColor','c')
plot(getnormal(Enzm),'c',ampl)
% plot(Enzmvy,'EdgeColor','c')
% plot(Enzmvxy,'EdgeColor','c')

hold off
grid on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')

figure('Name','Ellipsoid')
surf(E,'EdgeColor','k')
hold on
plot(getnormal(E),'k',ampl)
% surf(Evy,'EdgeColor','k')
% surf(Evxy,'EdgeColor','k')

surf(Enx,'EdgeColor','b')
plot(getnormal(Enx),'b',ampl)
% surf(Enxvy,'EdgeColor','b')
% surf(Enxvxy,'EdgeColor','b')

surf(Eny,'EdgeColor','r')
plot(getnormal(Eny),'r',ampl)
% surf(Enyvy,'EdgeColor','r')
% surf(Enyvxy,'EdgeColor','r')

surf(Enxyz,'EdgeColor','g')
plot(getnormal(Enxyz),'g',ampl)
% surf(Enxyzvy,'EdgeColor','g')
% surf(Enxyzvxy,'EdgeColor','g')

surf(Enzm,'EdgeColor','c')
plot(getnormal(Enzm),'c',ampl)
% surf(Enzmvy,'EdgeColor','c')
% surf(Enzmvxy,'EdgeColor','c')

hold off
grid on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')