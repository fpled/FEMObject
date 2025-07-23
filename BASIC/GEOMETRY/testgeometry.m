center = [0,0,0];

vx = [1,0];
vy = [0,1];
vxy = [1,1];

nx = [1,0,0];
ny = [0,1,0];
nz = [0,0,1];
nxy = [1,1,0];
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

Cnxy = CIRCLE(center(1),center(2),center(3),r,nxy(1),nxy(2),nxy(3));
Cnxyvy = CIRCLE(center(1),center(2),center(3),r,nxy(1),nxy(2),nxy(3),vy(1),vy(2));
Cnxyvxy = CIRCLE(center(1),center(2),center(3),r,nxy(1),nxy(2),nxy(3),vxy(1),vxy(2));

Cnxyz = CIRCLE(center(1),center(2),center(3),r,nxyz(1),nxyz(2),nxyz(3));
Cnxyzvy = CIRCLE(center(1),center(2),center(3),r,nxyz(1),nxyz(2),nxyz(3),vy(1),vy(2));
Cnxyzvxy = CIRCLE(center(1),center(2),center(3),r,nxyz(1),nxyz(2),nxyz(3),vxy(1),vxy(2));

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

plot(Cnxy,'EdgeColor','g')
plot(getnormal(Cnxy),'g',ampl)
% plot(Cnxyvy,'EdgeColor','g')
% plot(Cnxyvxy,'EdgeColor','g')

plot(Cnxyz,'EdgeColor','m')
plot(getnormal(Cnxyz),'m',ampl)
% plot(Cnxyzvy,'EdgeColor','m')
% plot(Cnxyzvxy,'EdgeColor','m')

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

Evy  = ELLIPSE(center(1),center(2),center(3),a,b,vy(1),vy(2));
Evxy = ELLIPSE(center(1),center(2),center(3),a,b,vxy(1),vxy(2));

Enx    = ELLIPSE(center(1),center(2),center(3),a,b,nx(1),nx(2),nx(3));
Enxvy  = ELLIPSE(center(1),center(2),center(3),a,b,nx(1),nx(2),nx(3),vy(1),vy(2));
Enxvxy = ELLIPSE(center(1),center(2),center(3),a,b,nx(1),nx(2),nx(3),vxy(1),vxy(2));

Eny    = ELLIPSE(center(1),center(2),center(3),a,b,ny(1),ny(2),ny(3));
Enyvy  = ELLIPSE(center(1),center(2),center(3),a,b,ny(1),ny(2),ny(3),vy(1),vy(2));
Enyvxy = ELLIPSE(center(1),center(2),center(3),a,b,ny(1),ny(2),ny(3),vxy(1),vxy(2));

Enxy    = ELLIPSE(center(1),center(2),center(3),a,b,nxy(1),nxy(2),nxy(3));
Enxyvy  = ELLIPSE(center(1),center(2),center(3),a,b,nxy(1),nxy(2),nxy(3),vy(1),vy(2));
Enxyvxy = ELLIPSE(center(1),center(2),center(3),a,b,nxy(1),nxy(2),nxy(3),vxy(1),vxy(2));

Enxyz    = ELLIPSE(center(1),center(2),center(3),a,b,nxyz(1),nxyz(2),nxyz(3));
Enxyzvy  = ELLIPSE(center(1),center(2),center(3),a,b,nxyz(1),nxyz(2),nxyz(3),vy(1),vy(2));
Enxyzvxy = ELLIPSE(center(1),center(2),center(3),a,b,nxyz(1),nxyz(2),nxyz(3),vxy(1),vxy(2));

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

plot(Enxy,'EdgeColor','g')
plot(getnormal(Enxy),'g',ampl)
% plot(Enxyvy,'EdgeColor','g')
% plot(Enxyvxy,'EdgeColor','g')

plot(Enxyz,'EdgeColor','m')
plot(getnormal(Enxyz),'m',ampl)
% plot(Enxyzvy,'EdgeColor','m')
% plot(Enxyzvxy,'EdgeColor','m')

hold off
grid on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')

%% CYLINDER
r = 0.3;
h = 1;
angle = 2*pi;
% angle = 4*pi/3;
ampl_plot = h/2;
ampl_surf = h;

CY = CYLINDER(center(1),center(2),center(3),r,h,angle);

CYvy  = CYLINDER(center(1),center(2),center(3),r,h,nz(1),nz(2),nz(3),vy(1),vy(2),angle);
CYvxy = CYLINDER(center(1),center(2),center(3),r,h,nz(1),nz(2),nz(3),vxy(1),vxy(2),angle);

CYnx    = CYLINDER(center(1),center(2),center(3),r,h,nx(1),nx(2),nx(3),angle);
CYnxvy  = CYLINDER(center(1),center(2),center(3),r,h,nx(1),nx(2),nx(3),vy(1),vy(2),angle);
CYnxvxy = CYLINDER(center(1),center(2),center(3),r,h,nx(1),nx(2),nx(3),vxy(1),vxy(2),angle);

CYny    = CYLINDER(center(1),center(2),center(3),r,h,ny(1),ny(2),ny(3),angle);
CYnyvy  = CYLINDER(center(1),center(2),center(3),r,h,ny(1),ny(2),ny(3),vy(1),vy(2),angle);
CYnyvxy = CYLINDER(center(1),center(2),center(3),r,h,ny(1),ny(2),ny(3),vxy(1),vxy(2),angle);

CYnxy    = CYLINDER(center(1),center(2),center(3),r,h,nxy(1),nxy(2),nxy(3),angle);
CYnxyvy  = CYLINDER(center(1),center(2),center(3),r,h,nxy(1),nxy(2),nxy(3),vy(1),vy(2),angle);
CYnxyvxy = CYLINDER(center(1),center(2),center(3),r,h,nxy(1),nxy(2),nxy(3),vxy(1),vxy(2),angle);

CYnxyz    = CYLINDER(center(1),center(2),center(3),r,h,nxyz(1),nxyz(2),nxyz(3),angle);
CYnxyzvy  = CYLINDER(center(1),center(2),center(3),r,h,nxyz(1),nxyz(2),nxyz(3),vy(1),vy(2),angle);
CYnxyzvxy = CYLINDER(center(1),center(2),center(3),r,h,nxyz(1),nxyz(2),nxyz(3),vxy(1),vxy(2),angle);

figure('Name','Cylinder')
plot(CY,'EdgeColor','k')
hold on
plot(getnormal(CY),'k',ampl_plot)
% plot(CYvy,'EdgeColor','k')
% plot(CYvxy,'EdgeColor','k')

plot(CYnx,'EdgeColor','b')
plot(getnormal(CYnx),'b',ampl_plot)
% plot(CYnxvy,'EdgeColor','b')
% plot(CYnxvxy,'EdgeColor','b')

plot(CYny,'EdgeColor','r')
plot(getnormal(CYny),'r',ampl_plot)
% plot(CYnyvy,'EdgeColor','r')
% plot(CYnyvxy,'EdgeColor','r')

plot(CYnxy,'EdgeColor','g')
plot(getnormal(CYnxy),'g',ampl_plot)
% plot(CYnxyvy,'EdgeColor','g')
% plot(CYnxyvxy,'EdgeColor','g')

plot(CYnxyz,'EdgeColor','m')
plot(getnormal(CYnxyz),'m',ampl_plot)
% plot(CYnxyzvy,'EdgeColor','m')
% plot(CYnxyzvxy,'EdgeColor','m')

hold off
grid on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')

figure('Name','Cylinder')
surf(CY,'EdgeColor','k')
hold on
plot(getnormal(CY),'k',ampl_surf)
% surf(CYvy,'EdgeColor','k')
% surf(CYvxy,'EdgeColor','k')

surf(CYnx,'EdgeColor','b')
plot(getnormal(CYnx),'b',ampl_surf)
% surf(CYnxvy,'EdgeColor','b')
% surf(CYnxvxy,'EdgeColor','b')

surf(CYny,'EdgeColor','r')
plot(getnormal(CYny),'r',ampl_surf)
% surf(CYnyvy,'EdgeColor','r')
% surf(CYnyvxy,'EdgeColor','r')

surf(CYnxy,'EdgeColor','g')
plot(getnormal(CYnxy),'g',ampl_surf)
% surf(CYnxyvy,'EdgeColor','g')
% surf(CYnxyvxy,'EdgeColor','g')

surf(CYnxyz,'EdgeColor','m')
plot(getnormal(CYnxyz),'m',ampl)
% surf(CYnxyzvy,'EdgeColor','m')
% surf(CYnxyzvxy,'EdgeColor','m')

hold off
grid on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')

%% SPHERE
r = 0.3;
ampl = r/2;

S = SPHERE(center(1),center(2),center(3),r);

Svy  = SPHERE(center(1),center(2),center(3),r,vy(1),vy(2));
Svxy = SPHERE(center(1),center(2),center(3),r,vxy(1),vxy(2));

Snx    = SPHERE(center(1),center(2),center(3),r,nx(1),nx(2),nx(3));
Snxvy  = SPHERE(center(1),center(2),center(3),r,nx(1),nx(2),nx(3),vy(1),vy(2));
Snxvxy = SPHERE(center(1),center(2),center(3),r,nx(1),nx(2),nx(3),vxy(1),vxy(2));

Sny    = SPHERE(center(1),center(2),center(3),r,ny(1),ny(2),ny(3));
Snyvy  = SPHERE(center(1),center(2),center(3),r,ny(1),ny(2),ny(3),vy(1),vy(2));
Snyvxy = SPHERE(center(1),center(2),center(3),r,ny(1),ny(2),ny(3),vxy(1),vxy(2));

Snxy    = SPHERE(center(1),center(2),center(3),r,nxy(1),nxy(2),nxy(3));
Snxyvy  = SPHERE(center(1),center(2),center(3),r,nxy(1),nxy(2),nxy(3),vy(1),vy(2));
Snxyvxy = SPHERE(center(1),center(2),center(3),r,nxy(1),nxy(2),nxy(3),vxy(1),vxy(2));

Snxyz    = SPHERE(center(1),center(2),center(3),r,nxyz(1),nxyz(2),nxyz(3));
Snxyzvy  = SPHERE(center(1),center(2),center(3),r,nxyz(1),nxyz(2),nxyz(3),vy(1),vy(2));
Snxyzvxy = SPHERE(center(1),center(2),center(3),r,nxyz(1),nxyz(2),nxyz(3),vxy(1),vxy(2));

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

plot(Snxy,'EdgeColor','g')
plot(getnormal(Snxy),'g',ampl)
% plot(Snxyvy,'EdgeColor','g')
% plot(Snxyvxy,'EdgeColor','g')

plot(Snxyz,'EdgeColor','m')
plot(getnormal(Snxyz),'m',ampl)
% plot(Snxyzvy,'EdgeColor','m')
% plot(Snxyzvxy,'EdgeColor','m')

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

surf(Snxy,'EdgeColor','g')
plot(getnormal(Snxy),'g',ampl)
% surf(Snxyvy,'EdgeColor','g')
% surf(Snxyvxy,'EdgeColor','g')

surf(Snxyz,'EdgeColor','m')
plot(getnormal(Snxyz),'m',ampl)
% surf(Snxyzvy,'EdgeColor','m')
% surf(Snxyzvxy,'EdgeColor','m')

hold off
grid on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')

%% ELLIPSOID
a = 0.5;
b = 0.2;
c = 0.1;
ampl = a/2;

EL = ELLIPSOID(center(1),center(2),center(3),a,b,c);

ELvy  = ELLIPSOID(center(1),center(2),center(3),a,b,c,vy(1),vy(2));
ELvxy = ELLIPSOID(center(1),center(2),center(3),a,b,c,vxy(1),vxy(2));

ELnx    = ELLIPSOID(center(1),center(2),center(3),a,b,c,nx(1),nx(2),nx(3));
ELnxvy  = ELLIPSOID(center(1),center(2),center(3),a,b,c,nx(1),nx(2),nx(3),vy(1),vy(2));
ELnxvxy = ELLIPSOID(center(1),center(2),center(3),a,b,c,nx(1),nx(2),nx(3),vxy(1),vxy(2));

ELny    = ELLIPSOID(center(1),center(2),center(3),a,b,c,ny(1),ny(2),ny(3));
ELnyvy  = ELLIPSOID(center(1),center(2),center(3),a,b,c,ny(1),ny(2),ny(3),vy(1),vy(2));
ELnyvxy = ELLIPSOID(center(1),center(2),center(3),a,b,c,ny(1),ny(2),ny(3),vxy(1),vxy(2));

ELnxy    = ELLIPSOID(center(1),center(2),center(3),a,b,c,nxy(1),nxy(2),nxy(3));
ELnxyvy  = ELLIPSOID(center(1),center(2),center(3),a,b,c,nxy(1),nxy(2),nxy(3),vy(1),vy(2));
ELnxyvxy = ELLIPSOID(center(1),center(2),center(3),a,b,c,nxy(1),nxy(2),nxy(3),vxy(1),vxy(2));

ELnxyz    = ELLIPSOID(center(1),center(2),center(3),a,b,c,nxyz(1),nxyz(2),nxyz(3));
ELnxyzvy  = ELLIPSOID(center(1),center(2),center(3),a,b,c,nxyz(1),nxyz(2),nxyz(3),vy(1),vy(2));
ELnxyzvxy = ELLIPSOID(center(1),center(2),center(3),a,b,c,nxyz(1),nxyz(2),nxyz(3),vxy(1),vxy(2));

figure('Name','Ellipsoid')
plot(EL,'EdgeColor','k')
hold on
plot(getnormal(EL),'k',ampl)
% plot(ELvy,'EdgeColor','k')
% plot(ELvxy,'EdgeColor','k')

plot(ELnx,'EdgeColor','b')
plot(getnormal(ELnx),'b',ampl)
% plot(ELnxvy,'EdgeColor','b')
% plot(ELnxvxy,'EdgeColor','b')

plot(ELny,'EdgeColor','r')
plot(getnormal(ELny),'r',ampl)
% plot(ELnyvy,'EdgeColor','r')
% plot(ELnyvxy,'EdgeColor','r')

plot(ELnxy,'EdgeColor','g')
plot(getnormal(ELnxy),'g',ampl)
% plot(ELnxyvy,'EdgeColor','g')
% plot(ELnxyvxy,'EdgeColor','g')

plot(ELnxyz,'EdgeColor','m')
plot(getnormal(ELnxyz),'m',ampl)
% plot(ELnxyzvy,'EdgeColor','m')
% plot(ELnxyzvxy,'EdgeColor','m')

hold off
grid on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')

figure('Name','Ellipsoid')
surf(EL,'EdgeColor','k')
hold on
plot(getnormal(EL),'k',ampl)
% surf(ELvy,'EdgeColor','k')
% surf(ELvxy,'EdgeColor','k')

surf(ELnx,'EdgeColor','b')
plot(getnormal(ELnx),'b',ampl)
% surf(ELnxvy,'EdgeColor','b')
% surf(ELnxvxy,'EdgeColor','b')

surf(ELny,'EdgeColor','r')
plot(getnormal(ELny),'r',ampl)
% surf(ELnyvy,'EdgeColor','r')
% surf(ELnyvxy,'EdgeColor','r')

surf(ELnxy,'EdgeColor','g')
plot(getnormal(ELnxy),'g',ampl)
% surf(ELnxyvy,'EdgeColor','g')
% surf(ELnxyvxy,'EdgeColor','g')

surf(ELnxyz,'EdgeColor','m')
plot(getnormal(ELnxyz),'m',ampl)
% surf(ELnxyzvy,'EdgeColor','m')
% surf(ELnxyzvxy,'EdgeColor','m')

hold off
grid on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')

%% TORUS
r1 = 0.5;
r2 = 0.1;
angle1 = pi/4;
angle2 = pi/2;
angle3 = pi;

T = TORUS(center(1),center(2),center(3),r1,r2);
T1 = TORUS(center(1),center(2),center(3),r1,r2,angle1);
T2 = TORUS(center(1),center(2),center(3),r1,r2,angle2);
T3 = TORUS(center(1),center(2),center(3),r1,r2,angle3);

figure('Name','Torus')
plot(T,'EdgeColor','k')
hold on
plot(T1,'EdgeColor','b')
plot(T2,'EdgeColor','r')
plot(T3,'EdgeColor','g')
hold off
grid on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')

figure('Name','Torus')
surf(T,'EdgeColor','k')
hold on
surf(T1,'EdgeColor','b')
surf(T2,'EdgeColor','r')
surf(T3,'EdgeColor','g')
hold off
grid on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')