%% Dimension 3
D = DOMAIN(3);
P1 = POINT([0.7,0.6,0.7]);
P2 = POINT([0.7,0.6,0.05]);

r = 0.3;
C1 = CIRCLE(0.5,0.5,0.0,r);
S1 = SPHERE(C1);
% S1 = SPHERE(0.5,0.5,0.0,r);
C2 = CIRCLE(0.5,0.5,0.5,r);
S2 = SPHERE(C2);
% S2 = SPHERE(0.5,0.5,0.5,r);

h1 = 1.0;
angle1 = 2*pi;
angle2 = 4*pi/3;
CY1 = CYLINDER(C1,h1,angle1);
% CY1 = CYLINDER(0.5,0.5,0.0,r,h1,angle1);
CY2 = CYLINDER(C1,h1,angle2);
% CY2 = CYLINDER(0.5,0.5,0.0,r,h1,angle2);

h2 = 0.5;
C3 = CIRCLE(0.5,0.5,0.25,r);
CY3 = CYLINDER(C3,h2,angle1);
% CY3 = CYLINDER(0.5,0.5,0.25,r,h2,angle1);
CY4 = CYLINDER(C3,h2,angle2);
% CY4 = CYLINDER(0.5,0.5,0.25,r,h2,angle2);

a = 0.4;
b = 0.2;
c = 0.1;
E = ELLIPSE(0.5,0.5,0.0,a,b);
EL1 = ELLIPSOID(E,c);
% EL1 = ELLIPSOID(0.5,0.5,0.0,a,b,c);
E2 = ELLIPSE(0.5,0.5,0.5,a,b);
EL2 = ELLIPSOID(E2,c);
% EL2 = ELLIPSOID(0.5,0.5,0.5,a,b,c);

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

%% Domain with cylindrical hole
fprintf('\n3D Domain with cylindrical hole\n');
fprintf('\n');
St = gmshDomainWithHole(D,CY3,0.1,0.02,'gmsh_domain_with_cylindrical_hole_tet');
% Sq = gmshDomainWithHole(D,CY3,0.1,0.02,'gmsh_domain_with_cylindrical_hole_cub',3,'recombine');

figure('Name','Domain with cylindrical hole')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

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

%% Domain with cylindrical inclusion
fprintf('\n3D Domain with cylindrical inclusion\n');
fprintf('\n');
St = gmshDomainWithInclusion(D,CY3,0.1,0.02,'gmsh_domain_with_cylindrical_inclusion_tet');
% Sq = gmshDomainWithInclusion(D,CY3,0.1,0.02,'gmsh_domain_with_cylindrical_inclusion_cub',3,'recombine');

figure('Name','Domain with cylindrical inclusion')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

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

%% Domain with single edge crack
fprintf('\n3D Domain with single edge crack\n');
St = gmshDomainWithSingleEdgeCrack(D,Q,0.1,0.02,'gmsh_domain_with_single_edge_crack_tet');
% Sq = gmshDomainWithSingleEdgeCrack(D,Q,0.1,0.02,'gmsh_domain_with_single_edge_crack_cub',3,'recombine');
Strefcrack = gmshDomainWithSingleEdgeCrack(D,Q,0.1,0.02,'gmsh_domain_with_single_edge_crack_refined_tet',3,'refinecrack');
% Sqrefcrack = gmshDomainWithSingleEdgeCrack(D,Q,0.1,0.02,'gmsh_domain_with_single_edge_crack_refined_cub',3,'recombine','refinecrack');

figure('Name','Domain with single edge crack')
clf
subplot(2,2,1)
plotparamelem(St,'group')
% subplot(2,2,2)
% plotparamelem(Sq,'group')
subplot(2,2,3)
plotparamelem(Strefcrack,'group')
% subplot(2,2,4)
% plotparamelem(Sqrefcrack,'group')

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
subplot(2,2,1)
plotparamelem(St,'group')
% subplot(2,2,2)
% plotparamelem(Sq,'group')
subplot(2,2,3)
plotparamelem(Strefcrack,'group')
% subplot(2,2,4)
% plotparamelem(Sqrefcrack,'group')

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
subplot(2,2,1)
plotparamelem(St,'group')
% subplot(2,2,2)
% plotparamelem(Sq,'group')
subplot(2,2,3)
plotparamelem(Strefcrack,'group')
% subplot(2,2,4)
% plotparamelem(Sqrefcrack,'group')

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
subplot(2,2,1)
plotparamelem(St,'group')
% subplot(2,2,2)
% plotparamelem(Sq,'group')
subplot(2,2,3)
plotparamelem(Strefcrack,'group')
% subplot(2,2,4)
% plotparamelem(Sqrefcrack,'group')

%% Domain with double edge crack (two symmetric or asymmetric edge cracks)
fprintf('\n3D Domain with double edge crack (two asymmetric edge cracks)\n');
fprintf('\n');
St = gmshDomainWithDoubleEdgeCrack(D,Qa,Qb,0.1,0.02,'gmsh_domain_with_double_edge_crack_tet');
% Sq = gmshDomainWithDoubleEdgeCrack(D,Qa,Qb,0.1,0.02,'gmsh_domain_with_double_edge_crack_cub',3,'recombine');
Strefcrack = gmshDomainWithDoubleEdgeCrack(D,Qa,Qb,0.1,0.02,'gmsh_domain_with_double_edge_crack_refined_tet',3,'refinecrack');
% Sqrefcrack = gmshDomainWithDoubleEdgeCrack(D,Qa,Qb,0.1,0.02,'gmsh_domain_with_double_edge_crack_refined_cub',3,'recombine','refinecrack');

figure('Name','Domain with double edge crack (two asymmetric edge cracks)')
clf
subplot(2,2,1)
plotparamelem(St,'group')
% subplot(2,2,2)
% plotparamelem(Sq,'group')
subplot(2,2,3)
plotparamelem(Strefcrack,'group')
% subplot(2,2,4)
% plotparamelem(Sqrefcrack,'group')

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
subplot(2,2,1)
plotparamelem(St,'group')
% subplot(2,2,2)
% plotparamelem(Sq,'group')
subplot(2,2,3)
plotparamelem(Strefcrack,'group')
% subplot(2,2,4)
% plotparamelem(Sqrefcrack,'group')

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
subplot(2,2,1)
plotparamelem(St,'group')
% subplot(2,2,2)
% plotparamelem(Sq,'group')
subplot(2,2,3)
plotparamelem(Strefcrack,'group')
% subplot(2,2,4)
% plotparamelem(Sqrefcrack,'group')

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
subplot(2,2,1)
plotparamelem(St,'group')
% subplot(2,2,2)
% plotparamelem(Sq,'group')
subplot(2,2,3)
plotparamelem(Strefcrack,'group')
% subplot(2,2,4)
% plotparamelem(Sqrefcrack,'group')

%% L-shaped panel
fprintf('\n3D L-shaped panel\n');
fprintf('\n');
St = gmshLshape(0.1,'gmsh_L_shape_tet',3);
% Sq = gmshLshape(0.1,'gmsh_L_shape_cub',3,'recombine');

figure('Name','L-shaped panel')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% L-shaped panel
fprintf('\n3D L-shaped panel\n');
fprintf('\n');
a = 250e-3; % half-length
b = 30e-3; % distance of applied load from the right edge
e = 100e-3; % thickness
St = gmshLshapedPanel(a,b,e,20e-3,5e-3,'gmsh_L_shaped_panel_tet',3);
% Sq = gmshLshapedPanel(a,b,e,20e-3,5e-3,'gmsh_L_shaped_panel_cub',3,'recombine');

figure('Name','L-shaped panel')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')