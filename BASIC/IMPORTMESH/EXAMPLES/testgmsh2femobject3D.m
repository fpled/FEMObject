% clc
clearvars
close all

%% 3D geometric objects
D = DOMAIN(3);
P1 = POINT([0.7,0.6,0.7]);
P2 = POINT([0.7,0.6,0.05]);
P3 = POINT([0.75,0.6,0.05]);

L3 = LINE([0.7,0.2,0.0],[0.9,0.1,0.0]);

r = 0.3;
C1 = CIRCLE(0.5,0.5,0.0,r);
C2 = CIRCLE(0.5,0.5,0.5,r);
S1 = SPHERE(C1); % S1 = SPHERE(0.5,0.5,0.0,r);
S2 = SPHERE(C2); % S2 = SPHERE(0.5,0.5,0.5,r);

n = [1,1,1];
v = [1,1];
C1rot = CIRCLE(0.5,0.5,0.0,r,n(1),n(2),n(3),v(1),v(2));
C2rot = CIRCLE(0.5,0.5,0.5,r,n(1),n(2),n(3),v(1),v(2));
S1rot = SPHERE(C1rot); % S1rot = SPHERE(0.5,0.5,0.0,r,n(1),n(2),n(3),v(1),v(2));
S2rot = SPHERE(C2rot); % S2rot = SPHERE(0.5,0.5,0.5,r,n(1),n(2),n(3),v(1),v(2));

h1 = 1.0;
angle1 = '2*Pi'; % angle1 = 2*pi;
% angle2 = '2*Pi - Pi/4'; % angle2 = 2*pi - pi/4;
angle2 = '4*Pi/3'; % angle2 = 4*pi/3;
CY1 = CYLINDER(C1,h1,angle1); % CY1 = CYLINDER(0.5,0.5,0.0,r,h1,angle1);
CY2 = CYLINDER(C1,h1,angle2); % CY2 = CYLINDER(0.5,0.5,0.0,r,h1,angle2);
CY1rot = CYLINDER(C1rot,h1,angle1); % CY1rot = CYLINDER(0.5,0.5,0.0,r,hn(1),n(2),n(3),v(1),v(2),1,angle1);
CY2rot = CYLINDER(C1rot,h1,angle2); % CY2rot = CYLINDER(0.5,0.5,0.0,r,hn(1),n(2),n(3),v(1),v(2),1,angle2);

h2 = 0.5;
C3 = CIRCLE(0.5,0.5,0.25,r);
CY3 = CYLINDER(C3,h2,angle1); % CY3 = CYLINDER(0.5,0.5,0.25,r,h2,angle1);
CY4 = CYLINDER(C3,h2,angle2); % CY4 = CYLINDER(0.5,0.5,0.25,r,h2,angle2);

C3rot = CIRCLE(0.5,0.5,0.25,r,n(1),n(2),n(3),v(1),v(2));
CY3rot = CYLINDER(C3rot,h2,angle1); % CY3rot = CYLINDER(0.5,0.5,0.25,r,h2,n(1),n(2),n(3),v(1),v(2),angle1);
CY4rot = CYLINDER(C3rot,h2,angle2); % CY4rot = CYLINDER(0.5,0.5,0.25,r,h2,n(1),n(2),n(3),v(1),v(2),angle2);

a = 0.4;
b = 0.2;
c = 0.1;
E1 = ELLIPSE(0.5,0.5,0.0,a,b);
E2 = ELLIPSE(0.5,0.5,0.5,a,b);
EL1 = ELLIPSOID(E1,c); % EL1 = ELLIPSOID(0.5,0.5,0.0,a,b,c);
EL2 = ELLIPSOID(E2,c); % EL2 = ELLIPSOID(0.5,0.5,0.5,a,b,c);

E1rot = ELLIPSE(0.5,0.5,0.0,a,b,n(1),n(2),n(3),v(1),v(2));
E2rot = ELLIPSE(0.5,0.5,0.5,a,b,n(1),n(2),n(3),v(1),v(2));
EL1rot = ELLIPSOID(E1rot,c); % EL1rot = ELLIPSOID(0.5,0.5,0.0,a,b,c,n(1),n(2),n(3),v(1),v(2));
EL2rot = ELLIPSOID(E2rot,c); % EL2rot = ELLIPSOID(0.5,0.5,0.5,a,b,c,n(1),n(2),n(3),v(1),v(2));

r1 = r;
r2 = 0.1;
angle1 = '2*Pi'; % angle1 = 2*pi;
% angle2 = '2*Pi - Pi/4'; % angle2 = 2*pi - pi/4;
% angle2 = '2*Pi - Pi/2'; % angle2 = 2*pi - pi/2;
angle2 = '4*Pi/3'; % angle2 = 4*pi/3;
T1 = TORUS(C1,r2,angle1); % T1 = TORUS(0.5,0.5,0.0,r1,r2,angle1);
T2 = TORUS(C1,r2,angle2); % T2 = TORUS(0.5,0.5,0.0,r1,r2,angle2);
T1rot = TORUS(C1rot,r2,angle1); % T1 = TORUS(0.5,0.5,0.0,r1,r2,n(1),n(2),n(3),v(1),v(2),angle1);
T2rot = TORUS(C1rot,r2,angle2); % T2 = TORUS(0.5,0.5,0.0,r1,r2,n(1),n(2),n(3),v(1),v(2),angle2);

T3 = TORUS(C2,r2,angle1); % T3 = TORUS(0.5,0.5,0.5,r1,r2,angle1);
T4 = TORUS(C2,r2,angle2); % T4 = TORUS(0.5,0.5,0.5,r1,r2,angle2);
T3rot = TORUS(C2rot,r2,angle1); % T3 = TORUS(0.5,0.5,0.5,r1,r2,n(1),n(2),n(3),v(1),v(2),angle1);
T4rot = TORUS(C2rot,r2,angle2); % T4 = TORUS(0.5,0.5,0.5,r1,r2,n(1),n(2),n(3),v(1),v(2),angle2);

a = 0.3; % crack length
b = 0.5; % height of crack surface
e = 1; % thickness
Q = QUADRANGLE([0.0,b,0.0],[a,b,0.0],[a,b,e],[0.0,b,e]);

b = 1; % domain size
a = 0.2*b; % crack length
e = 0.025*b; % eccentricity
Qa = QUADRANGLE([0.0,b/2+e,0.0],[a,b/2+e,0.0],[a,b/2+e,b],[0.0,b/2+e,b]);
Qb = QUADRANGLE([b,b/2-e,0.0],[b-a,b/2-e,0.0],[b-a,b/2-e,b],[b,b/2-e,b]);

%% Domain
fprintf('\n3D Domain\n');
fprintf('\n');
St = gmsh(D,0.1,'filename','gmsh_domain_tet');
% Sq = gmsh(D,0.1,'filename','gmsh_domain_cub','recombine');

figure('Name','Domain')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Domain using OpenCASCADE
fprintf('\n3D Domain using OpenCASCADE\n');
fprintf('\n');
St = gmsh(D,0.1,'filename','gmsh_domain_tet_occ','occ');
% Sq = gmsh(D,0.1,'filename','gmsh_domain_cub_occ','occ','recombine');

figure('Name','Domain using OpenCASCADE')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Extruded Domain
fprintf('\n3D Extruded Domain\n');
fprintf('\n');
St = gmsh(D,0.1,'filename','gmsh_domain_extrude_tet','extrude');
Sq = gmsh(D,0.1,'filename','gmsh_domain_extrude_cub','extrude','recombine');

figure('Name','Extruded Domain')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Extruded Domain using OpenCASCADE
fprintf('\n3D Extruded Domain using OpenCASCADE\n');
fprintf('\n');
St = gmsh(D,0.1,'filename','gmsh_domain_extrude_tet_occ','occ','extrude');
Sq = gmsh(D,0.1,'filename','gmsh_domain_extrude_cub_occ','occ','extrude','recombine');

figure('Name','Extruded Domain using OpenCASCADE')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with point
fprintf('\n3D Domain with point\n');
St = gmsh(D,P1,0.1,0.02,'filename','gmsh_domain_with_point_tet');
% Sq = gmsh(D,P1,0.1,0.02,'filename','gmsh_domain_with_point_cub','recombine');

figure('Name','Domain with point')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Domain with point using OpenCASCADE
fprintf('\n3D Domain with point using OpenCASCADE\n');
St = gmsh(D,P1,0.1,0.02,'filename','gmsh_domain_with_point_tet_occ','occ');
% Sq = gmsh(D,P1,0.1,0.02,'filename','gmsh_domain_with_point_cub_occ','occ','recombine');

figure('Name','Domain with point using OpenCASCADE')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Extruded Domain with point
fprintf('\n3D Extruded Domain with point\n');
St = gmsh(D,P1,0.1,0.02,'filename','gmsh_domain_with_point_extrude_tet','extrude');
% Sq = gmsh(D,P1,0.1,0.02,'filename','gmsh_domain_with_point_extrude_cub','extrude','recombine');

figure('Name','Extruded Domain with point')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Extruded Domain with point using OpenCASCADE
fprintf('\n3D Extruded Domain with point using OpenCASCADE\n');
St = gmsh(D,P1,0.1,0.02,'filename','gmsh_domain_with_point_extrude_tet_occ','occ','extrude');
% Sq = gmsh(D,P1,0.1,0.02,'filename','gmsh_domain_with_point_extrude_cub_occ','occ','extrude','recombine');

figure('Name','Extruded Domain with point using OpenCASCADE')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Circle
fprintf('\n2D Circle in 3D space\n');
fprintf('\n');
St = gmsh(C1,0.05,'filename','gmsh_circle_3D_tri');
Sq = gmsh(C1,0.05,'filename','gmsh_circle_3D_quad','recombine');

figure('Name','Circle')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle using OpenCASCADE
fprintf('\n2D Circle in 3D space using OpenCASCADE\n');
fprintf('\n');
St = gmsh(C1,0.05,'filename','gmsh_circle_3D_tri_occ','occ');
Sq = gmsh(C1,0.05,'filename','gmsh_circle_3D_quad_occ','occ','recombine');

figure('Name','Circle using OpenCASCADE')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Ellipse
fprintf('\n2D Ellipse in 3D space\n');
fprintf('\n');
St = gmsh(E1,0.05,'filename','gmsh_ellipse_3D_tri');
Sq = gmsh(E1,0.05,'filename','gmsh_ellipse_3D_quad','recombine');

figure('Name','Ellipse')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Ellipse using OpenCASCADE
fprintf('\n2D Ellipse in 3D space using OpenCASCADE\n');
fprintf('\n');
St = gmsh(E1,0.05,'filename','gmsh_ellipse_3D_tri_occ','occ');
Sq = gmsh(E1,0.05,'filename','gmsh_ellipse_3D_quad_occ','occ','recombine');

figure('Name','Ellipse using OpenCASCADE')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Sphere
fprintf('\n3D Sphere\n');
fprintf('\n');
St = gmsh(S1,0.05,'filename','gmsh_sphere_tet');
% Sq = gmsh(S1,0.05,'filename','gmsh_sphere_cub','recombine');

figure('Name','Sphere')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Sphere using OpenCASCADE
fprintf('\n3D Sphere using OpenCASCADE\n');
fprintf('\n');
St = gmsh(S1,0.05,'filename','gmsh_sphere_tet_occ','occ');
% Sq = gmsh(S1,0.05,'filename','gmsh_sphere_cub_occ','occ','recombine');

figure('Name','Sphere using OpenCASCADE')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Sphere with point
fprintf('\n3D Sphere with point\n');
fprintf('\n');
St = gmsh(S1,P2,0.05,0.01,'filename','gmsh_sphere_with_point_tet');
% Sq = gmsh(S1,P2,0.05,0.01,'filename','gmsh_sphere_with_point_cub','recombine');

figure('Name','Sphere with point')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Sphere with point using OpenCASCADE
fprintf('\n3D Sphere with point using OpenCASCADE\n');
fprintf('\n');
St = gmsh(S1,P2,0.05,0.01,'filename','gmsh_sphere_with_point_tet_occ','occ');
% Sq = gmsh(S1,P2,0.05,0.01,'filename','gmsh_sphere_with_point_cub_occ','occ','recombine');

figure('Name','Sphere with point using OpenCASCADE')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Ellipsoid
fprintf('\n3D Ellipsoid\n');
fprintf('\n');
St = gmsh(EL1,0.02,'filename','gmsh_ellipsoid_tet');
% Sq = gmsh(EL1,0.02,'filename','gmsh_ellipsoid_cub','recombine');

figure('Name','Ellipsoid')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Ellipsoid with point
fprintf('\n3D Ellipsoid with point\n');
fprintf('\n');
St = gmsh(EL1,P2,0.02,0.004,'filename','gmsh_ellipsoid_with_point_tet');
% Sq = gmsh(EL1,P2,0.02,0.004,'filename','gmsh_ellipsoid_with_point_cub','recombine');

figure('Name','Ellipsoid with point')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Full Cylinder
fprintf('\n3D Full Cylinder\n');
fprintf('\n');
St = gmsh(CY1,0.05,'filename','gmsh_cylinder_full_tet');
% Sq = gmsh(CY1,0.05,'filename','gmsh_cylinder_full_cub','recombine');

figure('Name','Full Cylinder')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Full Cylinder using OpenCASCADE
fprintf('\n3D Full Cylinder using OpenCASCADE\n');
fprintf('\n');
St = gmsh(CY1,0.05,'filename','gmsh_cylinder_full_tet_occ','occ');
% Sq = gmsh(CY1,0.05,'filename','gmsh_cylinder_full_cub_occ','occ','recombine');

figure('Name','Full Cylinder using OpenCASCADE')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Extruded Full Cylinder
fprintf('\n3D Extruded Full Cylinder\n');
fprintf('\n');
St = gmsh(CY1,0.05,'filename','gmsh_cylinder_full_extrude_tet','extrude');
Sq = gmsh(CY1,0.05,'filename','gmsh_cylinder_full_extrude_cub','extrude','recombine');

figure('Name','Extruded Full Cylinder')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Extruded Full Cylinder using OpenCASCADE
fprintf('\n3D Extruded Full Cylinder using OpenCASCADE\n');
fprintf('\n');
St = gmsh(CY1,0.05,'filename','gmsh_cylinder_full_extrude_tet_occ','occ','extrude');
Sq = gmsh(CY1,0.05,'filename','gmsh_cylinder_full_extrude_cub_occ','occ','extrude','recombine');

figure('Name','Extruded Full Cylinder using OpenCASCADE')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Full Cylinder with point
fprintf('\n3D Full Cylinder with point\n');
fprintf('\n');
St = gmsh(CY1,P1,0.05,0.01,'filename','gmsh_cylinder_full_with_point_tet');
% Sq = gmsh(CY1,P1,0.05,0.01,'filename','gmsh_cylinder_full_with_point_cub','recombine');

figure('Name','Full Cylinder with point')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Full Cylinder with point using OpenCASCADE
fprintf('\n3D Full Cylinder with point using OpenCASCADE\n');
fprintf('\n');
St = gmsh(CY1,P1,0.05,0.01,'filename','gmsh_cylinder_full_with_point_tet_occ','occ');
% Sq = gmsh(CY1,P1,0.05,0.01,'filename','gmsh_cylinder_full_with_point_cub_occ','occ','recombine');

figure('Name','Full Cylinder with point using OpenCASCADE')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Extruded Full Cylinder with point
fprintf('\n3D Extruded Full Cylinder with point\n');
fprintf('\n');
St = gmsh(CY1,P1,0.05,0.01,'filename','gmsh_cylinder_full_with_point_extrude_tet','extrude');
% Sq = gmsh(CY1,P1,0.05,0.01,'filename','gmsh_cylinder_full_with_point_extrude_cub','extrude','recombine');

figure('Name','Extruded Full Cylinder with point')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Extruded Full Cylinder with point using OpenCASCADE
fprintf('\n3D Extruded Full Cylinder with point using OpenCASCADE\n');
fprintf('\n');
St = gmsh(CY1,P1,0.05,0.01,'filename','gmsh_cylinder_full_with_point_extrude_tet_occ','occ','extrude');
% Sq = gmsh(CY1,P1,0.05,0.01,'filename','gmsh_cylinder_full_with_point_extrude_cub_occ','occ','recombine','extrude');

figure('Name','Extruded Full Cylinder with point using OpenCASCADE')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Partial Cylinder
fprintf('\n3D Partial Cylinder\n');
fprintf('\n');
St = gmsh(CY2,0.05,'filename','gmsh_cylinder_partial_tet');
% Sq = gmsh(CY2,0.05,'filename','gmsh_cylinder_partial_cub','recombine');

figure('Name','Partial Cylinder')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Partial Cylinder using OpenCASCADE
fprintf('\n3D Partial Cylinder using OpenCASCADE\n');
fprintf('\n');
St = gmsh(CY2,0.05,'filename','gmsh_cylinder_partial_tet_occ','occ');
% Sq = gmsh(CY2,0.05,'filename','gmsh_cylinder_partial_cub_occ','occ','recombine');

figure('Name','Partial Cylinder using OpenCASCADE')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Extruded Partial Cylinder
fprintf('\n3D Extruded Partial Cylinder\n');
fprintf('\n');
St = gmsh(CY2,0.05,'filename','gmsh_cylinder_partial_extrude_tet','extrude');
Sq = gmsh(CY2,0.05,'filename','gmsh_cylinder_partial_extrude_cub','extrude','recombine');

figure('Name','Extruded Partial Cylinder')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Extruded Partial Cylinder using OpenCASCADE
fprintf('\n3D Extruded Partial Cylinder using OpenCASCADE\n');
fprintf('\n');
St = gmsh(CY2,0.05,'filename','gmsh_cylinder_partial_extrude_tet_occ','occ','extrude');
% Sq = gmsh(CY2,0.05,'filename','gmsh_cylinder_partial_extrude_cub_occ','occ','extrude','recombine');

figure('Name','Extruded Partial Cylinder using OpenCASCADE')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Partial Cylinder with point
fprintf('\n3D Partial Cylinder with point\n');
fprintf('\n');
St = gmsh(CY2,P1,0.05,0.01,'filename','gmsh_cylinder_partial_with_point_tet');
% Sq = gmsh(CY2,P1,0.05,0.01,'filename','gmsh_cylinder_partial_with_point_cub','recombine');

figure('Name','Partial Cylinder with point')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Partial Cylinder with point using OpenCASCADE
fprintf('\n3D Partial Cylinder with point using OpenCASCADE\n');
fprintf('\n');
St = gmsh(CY2,P1,0.05,0.01,'filename','gmsh_cylinder_partial_with_point_tet_occ','occ');
% Sq = gmsh(CY2,P1,0.05,0.01,'filename','gmsh_cylinder_partial_with_point_cub_occ','occ','recombine');

figure('Name','Partial Cylinder with point using OpenCASCADE')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Extruded Partial Cylinder with point
fprintf('\n3D Extruded Partial Cylinder with point\n');
fprintf('\n');
St = gmsh(CY2,P1,0.05,0.01,'filename','gmsh_cylinder_partial_with_point_extrude_tet','extrude');
% Sq = gmsh(CY2,P1,0.05,0.01,'filename','gmsh_cylinder_partial_with_point_extrude_cub','extrude','recombine');

figure('Name','Extruded Partial Cylinder with point')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Extruded Partial Cylinder with point using OpenCASCADE
fprintf('\n3D Extruded Partial Cylinder with point using OpenCASCADE\n');
fprintf('\n');
St = gmsh(CY2,P1,0.05,0.01,'filename','gmsh_cylinder_partial_with_point_extrude_tet_occ','occ','extrude');
% Sq = gmsh(CY2,P1,0.05,0.01,'filename','gmsh_cylinder_partial_with_point_extrude_cub_occ','occ','recombine','extrude');

figure('Name','Extruded Partial Cylinder with point using OpenCASCADE')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Full Torus
fprintf('\n3D Full Torus\n');
fprintf('\n');
St = gmsh(T1,0.02,'filename','gmsh_torus_full_tet');
% Sq = gmsh(T1,0.02,'filename','gmsh_torus_full_cub','recombine');

figure('Name','Full Torus with point')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Full Torus using OpenCASCADE
fprintf('\n3D Full Torus using OpenCASCADE\n');
fprintf('\n');
St = gmsh(T1,0.02,'filename','gmsh_torus_full_tet_occ','occ');
% Sq = gmsh(T1,0.02,'filename','gmsh_torus_full_cub_occ','occ','recombine');

figure('Name','Full Torus using OpenCASCADE')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Extruded Full Torus
fprintf('\n3D Extruded Full Torus\n');
fprintf('\n');
St = gmsh(T1,0.02,'filename','gmsh_torus_full_extrude_tet','extrude');
Sq = gmsh(T1,0.02,'filename','gmsh_torus_full_extrude_cub','extrude','recombine');

figure('Name','Extruded Full Torus with point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Extruded Full Torus using OpenCASCADE
fprintf('\n3D Extruded Full Torus using OpenCASCADE\n');
fprintf('\n');
St = gmsh(T1,0.02,'filename','gmsh_torus_full_extrude_tet_occ','occ','extrude');
Sq = gmsh(T1,0.02,'filename','gmsh_torus_full_extrude_cub_occ','occ','extrude','recombine');

figure('Name','Extruded Full Torus using OpenCASCADE')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Full Torus with point
fprintf('\n3D Full Torus with point\n');
fprintf('\n');
St = gmsh(T1,P3,0.02,0.004,'filename','gmsh_torus_full_with_point_tet');
% Sq = gmsh(T1,P3,0.02,0.004,'filename','gmsh_torus_full_with_point_cub','recombine');

figure('Name','Full Torus with point')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Full Torus with point using OpenCASCADE
fprintf('\n3D Full Torus with point using OpenCASCADE\n');
fprintf('\n');
St = gmsh(T1,P3,0.02,0.004,'filename','gmsh_torus_full_with_point_tet_occ','occ');
% Sq = gmsh(T1,P3,0.02,0.004,'filename','gmsh_torus_full_with_point_cub_occ','occ','recombine');

figure('Name','Full Torus with point using OpenCASCADE')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Extruded Full Torus with point
fprintf('\n3D Extruded Full Torus with point\n');
fprintf('\n');
St = gmsh(T1,P3,0.02,0.004,'filename','gmsh_torus_full_with_point_extrude_tet','extrude');
% Sq = gmsh(T1,P3,0.02,0.004,'filename','gmsh_torus_full_with_point_extrude_cub','extrude','recombine');

figure('Name','Extruded Full Torus with point')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Extruded Full Torus with point using OpenCASCADE
fprintf('\n3D Extruded Full Torus with point using OpenCASCADE\n');
fprintf('\n');
St = gmsh(T1,P3,0.02,0.004,'filename','gmsh_torus_full_with_point_extrude_tet_occ','occ','extrude');
% Sq = gmsh(T1,P3,0.02,0.004,'filename','gmsh_torus_full_with_point_extrude_cub_occ','occ','extrude','recombine');

figure('Name','Extruded Full Torus with point using OpenCASCADE')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Partial Torus
fprintf('\n3D Partial Torus\n');
fprintf('\n');
St = gmsh(T2,0.02,'filename','gmsh_torus_partial_tet');
% Sq = gmsh(T2,0.02,'filename','gmsh_torus_partial_cub','recombine');

figure('Name','Partial Torus')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Partial Torus using OpenCASCADE
fprintf('\n3D Partial Torus using OpenCASCADE\n');
fprintf('\n');
St = gmsh(T2,0.02,'filename','gmsh_torus_partial_tet_occ','occ');
% Sq = gmsh(T2,0.02,'filename','gmsh_torus_partial_cub_occ','occ','recombine');

figure('Name','Partial Torus using OpenCASCADE')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Extruded Partial Torus
fprintf('\n3D Extruded Partial Torus\n');
fprintf('\n');
St = gmsh(T2,0.02,'filename','gmsh_torus_partial_extrude_tet','extrude');
Sq = gmsh(T2,0.02,'filename','gmsh_torus_partial_extrude_cub','extrude','recombine');

figure('Name','Extruded Partial Torus')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Extruded Partial Torus using OpenCASCADE
fprintf('\n3D Extruded Partial Torus using OpenCASCADE\n');
fprintf('\n');
St = gmsh(T2,0.02,'filename','gmsh_torus_partial_extrude_tet_occ','occ','extrude');
Sq = gmsh(T2,0.02,'filename','gmsh_torus_partial_extrude_cub_occ','occ','extrude','recombine');

figure('Name','Extruded Partial Torus using OpenCASCADE')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Partial Torus with point
fprintf('\n3D Partial Torus with point\n');
fprintf('\n');
St = gmsh(T2,P3,0.02,0.004,'filename','gmsh_torus_partial_with_point_tet');
% Sq = gmsh(T2,P3,0.02,0.004,'filename','gmsh_torus_partial_with_point_cub','recombine');

figure('Name','Partial Torus with point')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Partial Torus with point using OpenCASCADE
fprintf('\n3D Partial Torus with point using OpenCASCADE\n');
fprintf('\n');
St = gmsh(T2,P3,0.02,0.004,'filename','gmsh_torus_partial_with_point_tet_occ','occ');
% Sq = gmsh(T2,P3,0.02,0.004,'filename','gmsh_torus_partial_with_point_cub_occ','occ','recombine');

figure('Name','Partial Torus with point using OpenCASCADE')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Extruded Partial Torus with point
fprintf('\n3D Extruded Partial Torus with point\n');
fprintf('\n');
St = gmsh(T2,P3,0.02,0.004,'filename','gmsh_torus_partial_with_point_extrude_tet','extrude');
% Sq = gmsh(T2,P3,0.02,0.004,'filename','gmsh_torus_partial_with_point_extrude_cub','extrude','recombine');

figure('Name','Extruded Partial Torus with point')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Extruded Partial Torus with point using OpenCASCADE
fprintf('\n3D Extruded Partial Torus with point using OpenCASCADE\n');
fprintf('\n');
St = gmsh(T2,P3,0.02,0.004,'filename','gmsh_torus_partial_with_point_extrude_tet_occ','occ','extrude');
% Sq = gmsh(T2,P3,0.02,0.004,'filename','gmsh_torus_partial_with_point_extrude_cub_occ','occ','extrude','recombine');

figure('Name','Extruded Partial Torus with point using OpenCASCADE')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Extruded Domain with crack
fprintf('\n3D Extruded Domain with crack\n');
fprintf('\n');
St = gmshDomainWithHole(D,L3,0.1,0.02,'gmsh_domain_with_crack_extrude_tet',3,'extrude');
Sq = gmshDomainWithHole(D,L3,0.1,0.02,'gmsh_domain_with_crack_extrude_cub',3,'extrude','recombine');

figure('Name','Extruded Domain with crack line')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Extruded Domain with circular through-hole
fprintf('\n3D Domain with circular through-hole\n');
fprintf('\n');
St = gmshDomainWithHole(D,C1,0.1,0.02,'gmsh_domain_with_circular_hole_extrude_tet',3,'extrude');
Sq = gmshDomainWithHole(D,C1,0.1,0.02,'gmsh_domain_with_circular_hole_extrude_cub',3,'extrude','recombine');

figure('Name','Extruded Domain with circular through-hole')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with full cylindrical through-hole
fprintf('\n3D Extruded Domain with full cylindrical through-hole\n');
fprintf('\n');
St = gmshDomainWithHole(D,CY1,0.1,0.02,'gmsh_domain_with_full_cylindrical_hole_extrude_tet',3,'extrude');
Sq = gmshDomainWithHole(D,CY1,0.1,0.02,'gmsh_domain_with_full_cylindrical_hole_extrude_cub',3,'extrude','recombine');

figure('Name','Extruded Domain with full cylindrical through-hole')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with partial cylindrical hole
fprintf('\n3D Extruded Domain with partial cylindrical hole\n');
fprintf('\n');
St = gmshDomainWithHole(D,CY2,0.1,0.02,'gmsh_domain_with_partial_cylindrical_hole_extrude_tet',3,'extrude');
Sq = gmshDomainWithHole(D,CY2,0.1,0.02,'gmsh_domain_with_partial_cylindrical_hole_extrude_cub',3,'extrude','recombine');

figure('Name','Extruded Domain with partial cylindrical through-hole')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with spherical hole
fprintf('\n3D Domain with spherical hole\n');
fprintf('\n');
St = gmshDomainWithHole(D,S2,0.1,0.02,'gmsh_domain_with_spherical_hole_tet');
% Sq = gmshDomainWithHole(D,S2,0.1,0.02,'gmsh_domain_with_spherical_hole_cub',3,'recombine');

figure('Name','Domain with spherical hole')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Domain with ellipsoidal hole
fprintf('\n3D Domain with ellipsoidal hole\n');
fprintf('\n');
St = gmshDomainWithHole(D,EL2,0.1,0.02,'gmsh_domain_with_ellipsoidal_hole_tet');
% Sq = gmshDomainWithHole(D,EL2,0.1,0.02,'gmsh_domain_with_ellipsoidal_hole_cub',3,'recombine');

figure('Name','Domain with ellipsoidal hole')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Domain with full cylindrical hole
fprintf('\n3D Domain with full cylindrical hole\n');
fprintf('\n');
St = gmshDomainWithHole(D,CY3,0.1,0.02,'gmsh_domain_with_full_cylindrical_hole_tet');
% Sq = gmshDomainWithHole(D,CY3,0.1,0.02,'gmsh_domain_with_full_cylindrical_hole_cub',3,'recombine');

figure('Name','Domain with full cylindrical hole')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Domain with partial cylindrical hole
fprintf('\n3D Domain with partial cylindrical hole\n');
fprintf('\n');
St = gmshDomainWithHole(D,CY4,0.1,0.02,'gmsh_domain_with_partial_cylindrical_hole_tet');
% Sq = gmshDomainWithHole(D,CY4,0.1,0.02,'gmsh_domain_with_partial_cylindrical_hole_cub',3,'recombine');

figure('Name','Domain with partial cylindrical hole')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Domain with full toroidal hole
fprintf('\n3D Domain with full toroidal hole\n');
fprintf('\n');
St = gmshDomainWithHole(D,T3,0.1,0.02,'gmsh_domain_with_full_toroidal_hole_tet');
% Sq = gmshDomainWithHole(D,T3,0.1,0.02,'gmsh_domain_with_full_toroidal_hole_cub',3,'recombine');

figure('Name','Domain with full toroidal hole')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Domain with partial toroidal hole
fprintf('\n3D Domain with partial toroidal hole\n');
fprintf('\n');
St = gmshDomainWithHole(D,T4,0.1,0.02,'gmsh_domain_with_partial_toroidal_hole_tet');
% Sq = gmshDomainWithHole(D,T4,0.1,0.02,'gmsh_domain_with_partial_toroidal_hole_cub',3,'recombine');

figure('Name','Domain with partial toroidal hole')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Extruded Domain with line
fprintf('\n3D Extruded Domain with line\n');
fprintf('\n');
St = gmshDomainWithInclusion(D,L3,0.1,0.02,'gmsh_domain_with_line_extrude_tet',3,'extrude');
Sq = gmshDomainWithInclusion(D,L3,0.1,0.02,'gmsh_domain_with_line_extrude_cub',3,'extrude','recombine');

figure('Name','Extruded Domain with line')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Extruded Domain with circular through-inclusion
fprintf('\n3D Extruded Domain with circular through-inclusion\n');
fprintf('\n');
St = gmshDomainWithInclusion(D,C1,0.1,0.02,'gmsh_domain_with_circular_inclusion_extrude_tet',3,'extrude');
Sq = gmshDomainWithInclusion(D,C1,0.1,0.02,'gmsh_domain_with_circular_inclusion_extrude_cub',3,'extrude','recombine');

figure('Name','Extruded Domain with circular through-inclusion')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Extruded Domain with full cylindrical through-inclusion
fprintf('\n3D Extruded Domain with full cylindrical through-inclusion\n');
fprintf('\n');
St = gmshDomainWithInclusion(D,CY1,0.1,0.02,'gmsh_domain_with_full_cylindrical_inclusion_extrude_tet',3,'extrude');
Sq = gmshDomainWithInclusion(D,CY1,0.1,0.02,'gmsh_domain_with_full_cylindrical_inclusion_extrude_cub',3,'extrude','recombine');

figure('Name','Extruded Domain with full cylindrical through-inclusion')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Extruded Domain with partial cylindrical through-inclusion
fprintf('\n3D Extruded Domain with partial cylindrical through-inclusion\n');
fprintf('\n');
St = gmshDomainWithInclusion(D,CY2,0.1,0.02,'gmsh_domain_with_partial_cylindrical_inclusion_extrude_tet',3,'extrude');
Sq = gmshDomainWithInclusion(D,CY2,0.1,0.02,'gmsh_domain_with_partial_cylindrical_inclusion_extrude_cub',3,'extrude','recombine');

figure('Name','Extruded Domain with partial cylindrical through-inclusion')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with spherical inclusion
fprintf('\n3D Domain with spherical inclusion\n');
fprintf('\n');
St = gmshDomainWithInclusion(D,S2,0.1,0.02,'gmsh_domain_with_spherical_inclusion_tet');
% Sq = gmshDomainWithInclusion(D,S2,0.1,0.02,'gmsh_domain_with_spherical_inclusion_cub',3,'recombine');

figure('Name','Domain with spherical inclusion')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Domain with ellipsoidal inclusion
fprintf('\n3D Domain with ellipsoidal inclusion\n');
fprintf('\n');
St = gmshDomainWithInclusion(D,EL2,0.1,0.02,'gmsh_domain_with_ellipsoidal_inclusion_tet');
% Sq = gmshDomainWithInclusion(D,EL2,0.1,0.02,'gmsh_domain_with_ellipsoidal_inclusion_cub',3,'recombine');

figure('Name','Domain with ellipsoidal inclusion')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Domain with full cylindrical inclusion
fprintf('\n3D Domain with full cylindrical inclusion\n');
fprintf('\n');
St = gmshDomainWithInclusion(D,CY3,0.1,0.02,'gmsh_domain_with_full_cylindrical_inclusion_tet');
% Sq = gmshDomainWithInclusion(D,CY3,0.1,0.02,'gmsh_domain_with_full_cylindrical_inclusion_cub',3,'recombine');

figure('Name','Domain with full cylindrical inclusion')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Domain with partial cylindrical inclusion
fprintf('\n3D Domain with partial cylindrical inclusion\n');
fprintf('\n');
St = gmshDomainWithInclusion(D,CY4,0.1,0.02,'gmsh_domain_with_partial_cylindrical_inclusion_tet');
% Sq = gmshDomainWithInclusion(D,CY4,0.1,0.02,'gmsh_domain_with_partial_cylindrical_inclusion_cub',3,'recombine');

figure('Name','Domain with partial cylindrical inclusion')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Domain with full toroidal inclusion
fprintf('\n3D Domain with full toroidal inclusion\n');
fprintf('\n');
St = gmshDomainWithInclusion(D,T3,0.1,0.02,'gmsh_domain_with_full_toroidal_inclusion_tet');
% Sq = gmshDomainWithInclusion(D,T3,0.1,0.02,'gmsh_domain_with_full_toroidal_inclusion_cub',3,'recombine');

figure('Name','Domain with full toroidal inclusion')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Domain with partial toroidal inclusion
fprintf('\n3D Domain with partial toroidal inclusion\n');
fprintf('\n');
St = gmshDomainWithInclusion(D,T4,0.1,0.02,'gmsh_domain_with_partial_toroidal_inclusion_tet');
% Sq = gmshDomainWithInclusion(D,T4,0.1,0.02,'gmsh_domain_with_partial_toroidal_inclusion_cub',3,'recombine');

figure('Name','Domain with partial toroidal inclusion')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Domain with single edge crack
fprintf('\n3D Domain with single edge crack\n');
St = gmshDomainWithSingleEdgeCrack(D,Q,0.1,0.02,'gmsh_domain_with_single_edge_crack_tet');
% Sq = gmshDomainWithSingleEdgeCrack(D,Q,0.1,0.02,'gmsh_domain_with_single_edge_crack_cub',3,'recombine');
Strefcrack = gmshDomainWithSingleEdgeCrack(D,Q,0.1,0.02,'gmsh_domain_with_single_edge_crack_refined_tet',3,'refinecrack');
% Sqrefcrack = gmshDomainWithSingleEdgeCrack(D,Q,0.1,0.02,'gmsh_domain_with_single_edge_crack_refined_cub',3,'recombine','refinecrack');

figure('Name','Domain with single edge crack')
clf
% subplot(2,2,1)
% plotparamelem(St,'group')
% subplot(2,2,2)
% plotparamelem(Sq,'group')
% subplot(2,2,3)
% plotparamelem(Strefcrack,'group')
% subplot(2,2,4)
% plotparamelem(Sqrefcrack,'group')
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Strefcrack,'group')

%% Domain with single edge circular notch
fprintf('\n3D Domain with single circular notch\n');
fprintf('\n');
c = max(getsize(D))/50; % notch width
St = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.02,'gmsh_domain_with_single_edge_circular_notch_tet',3,'c');
% Sq = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.02,'gmsh_domain_with_single_edge_circular_notch_cub',3,'c','recombine');
Strefcrack = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.02,'gmsh_domain_with_single_edge_circular_notch_refined_tet',3,'c','refinecrack');
% Sqrefcrack = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.02,'gmsh_domain_with_single_edge_circular_notch_refined_cub',3,'c','recombine','refinecrack');

figure('Name','Domain with single edge circular notch')
clf
% subplot(2,2,1)
% plotparamelem(St,'group')
% subplot(2,2,2)
% plotparamelem(Sq,'group')
% subplot(2,2,3)
% plotparamelem(Strefcrack,'group')
% subplot(2,2,4)
% plotparamelem(Sqrefcrack,'group')
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Strefcrack,'group')

%% Domain with single edge rectangular notch
fprintf('\n3D Domain with single rectangular notch\n');
fprintf('\n');
c = max(getsize(D))/50; % notch width
St = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.02,'gmsh_domain_with_single_edge_rectangular_notch_tet',3,'r');
% Sq = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.02,'gmsh_domain_with_single_edge_rectangular_notch_cub',3,'r','recombine');
Strefcrack = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.02,'gmsh_domain_with_single_edge_rectangular_notch_refined_tet',3,'r','refinecrack');
% Sqrefcrack = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.02,'gmsh_domain_with_single_edge_rectangular_notch_refined_cub',3,'r','recombine','refinecrack');

figure('Name','Domain with single edge rectangular notch')
clf
% subplot(2,2,1)
% plotparamelem(St,'group')
% subplot(2,2,2)
% plotparamelem(Sq,'group')
% subplot(2,2,3)
% plotparamelem(Strefcrack,'group')
% subplot(2,2,4)
% plotparamelem(Sqrefcrack,'group')
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Strefcrack,'group')

%% Domain with single edge V notch
fprintf('\n3D Domain with single V notch\n');
fprintf('\n');
c = max(getsize(D))/50; % notch width
St = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.02,'gmsh_domain_with_single_edge_V_notch_tet',3,'v');
% Sq = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.02,'gmsh_domain_with_single_edge_V_notch_cub',3,'v','recombine');
Strefcrack = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.02,'gmsh_domain_with_single_edge_V_notch_refined_tet',3,'v','refinecrack');
% Sqrefcrack = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.02,'gmsh_domain_with_single_edge_V_notch_refined_cub',3,'v','recombine','refinecrack');

figure('Name','Domain with single edge V notch')
clf
% subplot(2,2,1)
% plotparamelem(St,'group')
% subplot(2,2,2)
% plotparamelem(Sq,'group')
% subplot(2,2,3)
% plotparamelem(Strefcrack,'group')
% subplot(2,2,4)
% plotparamelem(Sqrefcrack,'group')
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Strefcrack,'group')

%% Domain with double edge crack (two symmetric or asymmetric edge cracks)
fprintf('\n3D Domain with double edge crack (two asymmetric edge cracks)\n');
fprintf('\n');
St = gmshDomainWithDoubleEdgeCrack(D,Qa,Qb,0.1,0.02,'gmsh_domain_with_double_edge_crack_tet');
% Sq = gmshDomainWithDoubleEdgeCrack(D,Qa,Qb,0.1,0.02,'gmsh_domain_with_double_edge_crack_cub',3,'recombine');
Strefcrack = gmshDomainWithDoubleEdgeCrack(D,Qa,Qb,0.1,0.02,'gmsh_domain_with_double_edge_crack_refined_tet',3,'refinecrack');
% Sqrefcrack = gmshDomainWithDoubleEdgeCrack(D,Qa,Qb,0.1,0.02,'gmsh_domain_with_double_edge_crack_refined_cub',3,'recombine','refinecrack');

figure('Name','Domain with double edge crack (two asymmetric edge cracks)')
clf
% subplot(2,2,1)
% plotparamelem(St,'group')
% subplot(2,2,2)
% plotparamelem(Sq,'group')
% subplot(2,2,3)
% plotparamelem(Strefcrack,'group')
% subplot(2,2,4)
% plotparamelem(Sqrefcrack,'group')
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Strefcrack,'group')

%% Domain with double edge circular notch (two symmetric or asymmetric edge circular notches)
fprintf('\n3D Domain with double edge circular notch (two asymmetric edge circular notches)\n');
fprintf('\n');
c = max(getsize(D))/50; % notch width
St = gmshDomainWithDoubleEdgeNotch(D,Qa,Qb,c,0.1,0.02,'gmsh_domain_with_double_edge_circular_notch_tet',3,'c');
% Sq = gmshDomainWithDoubleEdgeNotch(D,Qa,Qb,c,0.1,0.02,'gmsh_domain_with_double_edge_circular_notch_cub',3,'c','recombine');
Strefcrack = gmshDomainWithDoubleEdgeNotch(D,Qa,Qb,c,0.1,0.02,'gmsh_domain_with_double_edge_circular_notch_refined_tet',3,'c','refinecrack');
% Sqrefcrack = gmshDomainWithDoubleEdgeNotch(D,Qa,Qb,c,0.1,0.02,'gmsh_domain_with_double_edge_circular_notch_refined_cub',3,'c','recombine','refinecrack');

figure('Name','Domain with double edge circular notch (two asymmetric edge circular notches)')
clf
% subplot(2,2,1)
% plotparamelem(St,'group')
% subplot(2,2,2)
% plotparamelem(Sq,'group')
% subplot(2,2,3)
% plotparamelem(Strefcrack,'group')
% subplot(2,2,4)
% plotparamelem(Sqrefcrack,'group')
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Strefcrack,'group')

%% Domain with double edge rectangular notch (two symmetric or asymmetric edge rectangular notches)
fprintf('\n3D Domain with double edge rectangular notch (two asymmetric edge rectangular notches)\n');
fprintf('\n');
c = max(getsize(D))/50; % notch width
St = gmshDomainWithDoubleEdgeNotch(D,Qa,Qb,c,0.1,0.02,'gmsh_domain_with_double_edge_rectangular_notch_tet',3,'r');
% Sq = gmshDomainWithDoubleEdgeNotch(D,Qa,Qb,c,0.1,0.02,'gmsh_domain_with_double_edge_rectangular_notch_cub',3,'r','recombine');
Strefcrack = gmshDomainWithDoubleEdgeNotch(D,Qa,Qb,c,0.1,0.02,'gmsh_domain_with_double_edge_rectangular_notch_refined_tet',3,'r','refinecrack');
% Sqrefcrack = gmshDomainWithDoubleEdgeNotch(D,Qa,Qb,c,0.1,0.02,'gmsh_domain_with_double_edge_rectangular_notch_refined_cub',3,'r','recombine','refinecrack');

figure('Name','Domain with double edge rectangular notch (two asymmetric edge rectangular notches)')
clf
% subplot(2,2,1)
% plotparamelem(St,'group')
% subplot(2,2,2)
% plotparamelem(Sq,'group')
% subplot(2,2,3)
% plotparamelem(Strefcrack,'group')
% subplot(2,2,4)
% plotparamelem(Sqrefcrack,'group')
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Strefcrack,'group')

%% Domain with double edge V notch (two symmetric or asymmetric edge V notches)
fprintf('\n3D Domain with double edge V notch (two asymmetric edge V notches)\n');
fprintf('\n');
c = max(getsize(D))/50; % notch width
St = gmshDomainWithDoubleEdgeNotch(D,Qa,Qb,c,0.1,0.02,'gmsh_domain_with_double_edge_V_notch_tet',3,'v');
% Sq = gmshDomainWithDoubleEdgeNotch(D,Qa,Qb,c,0.1,0.02,'gmsh_domain_with_double_edge_V_notch_cub',3,'v','recombine');
Strefcrack = gmshDomainWithDoubleEdgeNotch(D,Qa,Qb,c,0.1,0.02,'gmsh_domain_with_double_edge_V_notch_refined_tet',3,'v','refinecrack');
% Sqrefcrack = gmshDomainWithDoubleEdgeNotch(D,Qa,Qb,c,0.1,0.02,'gmsh_domain_with_double_edge_V_notch_refined_cub',3,'v','recombine','refinecrack');

figure('Name','Domain with double edge V notch (two asymmetric edge V notches)')
clf
% subplot(2,2,1)
% plotparamelem(St,'group')
% subplot(2,2,2)
% plotparamelem(Sq,'group')
% subplot(2,2,3)
% plotparamelem(Strefcrack,'group')
% subplot(2,2,4)
% plotparamelem(Sqrefcrack,'group')
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Strefcrack,'group')

%% L-shaped panel
fprintf('\n3D L-shaped panel\n');
fprintf('\n');
St = gmshLshape(0.1,'gmsh_L_shape_tet',3);
Sq = gmshLshape(0.1,'gmsh_L_shape_cub',3,'recombine');

figure('Name','L-shaped panel')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% L-shaped panel
fprintf('\n3D L-shaped panel\n');
fprintf('\n');
a = 250e-3; % half-length
b = 30e-3; % distance of applied load from the right edge
e = 100e-3; % thickness
St = gmshLshapedPanel(a,b,e,20e-3,5e-3,'gmsh_L_shaped_panel_tet',3);
Sq = gmshLshapedPanel(a,b,e,20e-3,5e-3,'gmsh_L_shaped_panel_cub',3,'recombine');

figure('Name','L-shaped panel')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')