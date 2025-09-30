% clc
clearvars
close all

%% 2D geometric objects
D = DOMAIN(2);
P1 = POINT([0.3,0.5]);
P2 = POINT([0.1,0.2]);

r1 = 0.3;
r2 = 0.7;
C1 = CIRCLE(0.5,0.5,r1);
C2 = CIRCLE(0.5,0.5,r2);
% C1 = CIRCLE(0.5,0.5,r1,1,1);
% C2 = CIRCLE(0.5,0.5,r2,1,1);

L1 = LINE([0.0,0.5],P1);
L2 = LINE([0.2,0.3],[0.6,0.7]);
L3 = LINE([0.7,0.2],[0.9,0.1]);

a = 0.4;
b = 0.2;
E1 = ELLIPSE(0.5,0.5,a,b);
% E1 = ELLIPSE(0.5,0.5,a,b,1,1);

b = 1; % domain length
a = 0.2*b; % crack length
e = 0.025*b; % eccentricity
La = LINE([0.0,b/2+e],[a,b/2+e]);
Lb = LINE([b-a,b/2-e],[b,b/2-e]);

D1 = DOMAIN(2,[0.9,0.45],[1.0,0.55]);
D2 = DOMAIN(2,[0.51,0.45],[0.61,0.55]);
D3 = DOMAIN(2,[0.1,0.45],[0.2,0.55]);

%% Domain
fprintf('\n2D Domain\n');
fprintf('\n');
St = gmsh(D,0.05,'filename','gmsh_domain_tri');
Sq = gmsh(D,0.05,'filename','gmsh_domain_quad','recombine');

figure('Name','Domain')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain using OpenCASCADE
fprintf('\n2D Domain using OpenCASCADE\n');
fprintf('\n');
St = gmsh(D,0.05,'filename','gmsh_domain_tri_occ','occ');
Sq = gmsh(D,0.05,'filename','gmsh_domain_quad_occ','occ','recombine');

figure('Name','Domain using OpenCASCADE')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with point
fprintf('\n2D Domain with point\n');
fprintf('\n');
St = gmsh(D,P1,0.05,0.01,'filename','gmsh_domain_with_point_tri');
Sq = gmsh(D,P1,0.05,0.01,'filename','gmsh_domain_with_point_quad','recombine');

figure('Name','Domain with point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with point using OpenCASCADE
fprintf('\n2D Domain with point using OpenCASCADE\n');
fprintf('\n');
St = gmsh(D,P1,0.05,0.01,'filename','gmsh_domain_with_point_tri_occ','occ');
Sq = gmsh(D,P1,0.05,0.01,'filename','gmsh_domain_with_point_quad_occ','occ','recombine');

figure('Name','Domain with point using OpenCASCADE')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle
fprintf('\n2D Circle\n');
fprintf('\n');
St = gmsh(C1,0.05,'filename','gmsh_circle_tri');
Sq = gmsh(C1,0.05,'filename','gmsh_circle_quad','recombine');

figure('Name','Circle')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle using OpenCASCADE
fprintf('\n2D Circle using OpenCASCADE\n');
fprintf('\n');
St = gmsh(C1,0.05,'filename','gmsh_circle_tri_occ','occ');
Sq = gmsh(C1,0.05,'filename','gmsh_circle_quad_occ','occ','recombine');

figure('Name','Circle using OpenCASCADE')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with point
fprintf('\n2D Circle with point\n');
fprintf('\n');
St = gmsh(C1,P1,0.05,0.005,'filename','gmsh_circle_with_point_tri');
Sq = gmsh(C1,P1,0.05,0.005,'filename','gmsh_circle_with_point_quad','recombine');

figure('Name','Circle with point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with point using OpenCASCADE
fprintf('\n2D Circle with point using OpenCASCADE\n');
fprintf('\n');
St = gmsh(C1,P1,0.05,0.005,'filename','gmsh_circle_with_point_tri_occ','occ');
Sq = gmsh(C1,P1,0.05,0.005,'filename','gmsh_circle_with_point_quad_occ','occ','recombine');

figure('Name','Circle with point using OpenCASCADE')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Ellipse
fprintf('\n2D Ellipse\n');
fprintf('\n');
St = gmsh(E1,0.05,'filename','gmsh_ellipse_tri');
Sq = gmsh(E1,0.05,'filename','gmsh_ellipse_quad','recombine');

figure('Name','Ellipse')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Ellipse using OpenCASCADE
fprintf('\n2D Ellipse using OpenCASCADE\n');
fprintf('\n');
St = gmsh(E1,0.05,'filename','gmsh_ellipse_tri_occ','occ');
Sq = gmsh(E1,0.05,'filename','gmsh_ellipse_quad_occ','occ','recombine');

figure('Name','Ellipse using OpenCASCADE')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Ellipse with point
fprintf('\n2D Ellipse with point\n');
fprintf('\n');
St = gmsh(E1,P1,0.05,0.005,'filename','gmsh_ellipse_with_point_tri');
Sq = gmsh(E1,P1,0.05,0.005,'filename','gmsh_ellipse_with_point_quad','recombine');

figure('Name','Ellipse with point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Ellipse with point using OpenCASCADE
fprintf('\n2D Ellipse with point using OpenCASCADE\n');
fprintf('\n');
St = gmsh(E1,P1,0.05,0.005,'filename','gmsh_ellipse_with_point_tri_occ','occ');
Sq = gmsh(E1,P1,0.05,0.005,'filename','gmsh_ellipse_with_point_quad_occ','occ','recombine');

figure('Name','Ellipse with point using OpenCASCADE')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with circular hole
fprintf('\n2D Domain with circular hole\n');
fprintf('\n');
St = gmshDomainWithHole(D,C1,0.05,0.02,'gmsh_domain_with_circular_hole_tri');
Sq = gmshDomainWithHole(D,C1,0.05,0.02,'gmsh_domain_with_circular_hole_quad',2,'recombine');

figure('Name','Domain with circular hole')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with ellipsoidal hole
fprintf('\n2D Domain with ellipsoidal hole\n');
fprintf('\n');
St = gmshDomainWithHole(D,E1,0.05,0.02,'gmsh_domain_with_ellipsoidal_hole_tri');
Sq = gmshDomainWithHole(D,E1,0.05,0.02,'gmsh_domain_with_ellipsoidal_hole_quad',2,'recombine');

figure('Name','Domain with ellipsoidal hole')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with circular hole and point
fprintf('\n2D Domain with circular hole and point\n');
fprintf('\n');
St = gmshDomainWithHole(D,{C1,P2},0.05,[0.02,0.005],'gmsh_domain_with_circular_hole_point_tri');
Sq = gmshDomainWithHole(D,{C1,P2},0.05,[0.02,0.005],'gmsh_domain_with_circular_hole_point_quad',2,'recombine');

figure('Name','Domain with circular hole and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with ellipsoidal hole and point
fprintf('\n2D Domain with ellipsoidal hole and point\n');
fprintf('\n');
St = gmshDomainWithHole(D,{E1,P2},0.05,[0.02,0.005],'gmsh_domain_with_ellipsoidal_hole_point_tri');
Sq = gmshDomainWithHole(D,{E1,P2},0.05,[0.02,0.005],'gmsh_domain_with_ellipsoidal_hole_point_quad',2,'recombine');

figure('Name','Domain with ellipsoidal hole and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with circular hole, crack and point
fprintf('\n2D Domain with circular hole, crack and point\n');
fprintf('\n');
St = gmshDomainWithHole(D,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_circular_hole_crack_point_tri');
Sq = gmshDomainWithHole(D,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_circular_hole_crack_point_quad',2,'recombine');

figure('Name','Domain with circular hole, crack and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with ellipsoidal hole, crack and point
fprintf('\n2D Domain with ellipsoidal hole, crack and point\n');
fprintf('\n');
St = gmshDomainWithHole(D,{E1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_ellipsoidal_hole_crack_point_tri');
Sq = gmshDomainWithHole(D,{E1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_ellipsoidal_hole_crack_point_quad',2,'recombine');

figure('Name','Domain with ellipsoidal hole, crack and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with circular hole, line and point
fprintf('\n2D Domain with circular hole, line and point\n');
fprintf('\n');
St = gmshDomainWithHole(D,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_circular_hole_line_point_tri',2,'noduplicate');
Sq = gmshDomainWithHole(D,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_circular_hole_line_point_quad',2,'noduplicate','recombine');

figure('Name','Domain with circular hole, line and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with ellipsoidal hole, line and point
fprintf('\n2D Domain with ellipsoidal hole, line and point\n');
fprintf('\n');
St = gmshDomainWithHole(D,{E1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_ellipsoidal_hole_line_point_tri',2,'noduplicate');
Sq = gmshDomainWithHole(D,{E1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_ellipsoidal_hole_line_point_quad',2,'noduplicate','recombine');

figure('Name','Domain with ellipsoidal hole, line and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with circular inclusion
fprintf('\n2D Domain with circular inclusion\n');
fprintf('\n');
St = gmshDomainWithInclusion(D,C1,0.05,0.02,'gmsh_domain_with_circular_inclusion_tri');
Sq = gmshDomainWithInclusion(D,C1,0.05,0.02,'gmsh_domain_with_circular_inclusion_quad',2,'recombine');

figure('Name','Domain with circular inclusion')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with ellipsoidal inclusion
fprintf('\n2D Domain with ellipsoidal inclusion\n');
fprintf('\n');
St = gmshDomainWithInclusion(D,E1,0.05,0.02,'gmsh_domain_with_ellipsoidal_inclusion_tri');
Sq = gmshDomainWithInclusion(D,E1,0.05,0.02,'gmsh_domain_with_ellipsoidal_inclusion_quad',2,'recombine');

figure('Name','Domain with ellipsoidal inclusion')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with circular inclusion and point
fprintf('\n2D Domain with circular inclusion and point\n');
fprintf('\n');
St = gmshDomainWithInclusion(D,{C1,P2},0.05,[0.02,0.005],'gmsh_domain_with_circular_inclusion_point_tri');
Sq = gmshDomainWithInclusion(D,{C1,P2},0.05,[0.02,0.005],'gmsh_domain_with_circular_inclusion_point_quad',2,'recombine');

figure('Name','Domain with circular inclusion and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with ellipsoidal inclusion and point
fprintf('\n2D Domain with ellipsoidal inclusion and point\n');
fprintf('\n');
St = gmshDomainWithInclusion(D,{E1,P2},0.05,[0.02,0.005],'gmsh_domain_with_ellipsoidal_inclusion_point_tri');
Sq = gmshDomainWithInclusion(D,{E1,P2},0.05,[0.02,0.005],'gmsh_domain_with_ellipsoidal_inclusion_point_quad',2,'recombine');

figure('Name','Domain with ellipsoidal inclusion and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with circular inclusion, line and point
fprintf('\n2D Domain with circular inclusion, line and point\n');
fprintf('\n');
St = gmshDomainWithInclusion(D,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_circular_inclusion_line_point_tri');
Sq = gmshDomainWithInclusion(D,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_circular_inclusion_line_point_quad',2,'recombine');

figure('Name','Domain with circular inclusion, line and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with ellipsoidal inclusion, line and point
fprintf('\n2D Domain with ellipsoidal inclusion, line and point\n');
fprintf('\n');
St = gmshDomainWithInclusion(D,{E1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_ellipsoidal_inclusion_line_point_tri');
Sq = gmshDomainWithInclusion(D,{E1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_ellipsoidal_inclusion_line_point_quad',2,'recombine');

figure('Name','Domain with ellipsoidal inclusion, line and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with circular hole
fprintf('\n2D Circle with circular hole\n');
fprintf('\n');
St = gmshCircleWithHole(C2,C1,0.05,0.02,'gmsh_circle_with_circular_hole_tri');
Sq = gmshCircleWithHole(C2,C1,0.05,0.02,'gmsh_circle_with_circular_hole_quad',2,'recombine');

figure('Name','Circle with circular hole')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with ellipsoidal hole
fprintf('\n2D Circle with ellipsoidal hole\n');
fprintf('\n');
St = gmshCircleWithHole(C2,E1,0.05,0.02,'gmsh_circle_with_ellipsoidal_hole_tri');
Sq = gmshCircleWithHole(C2,E1,0.05,0.02,'gmsh_circle_with_ellipsoidal_hole_quad',2,'recombine');

figure('Name','Circle with ellipsoidal hole')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with circular hole and point
fprintf('\n2D Circle with circular hole and point\n');
fprintf('\n');
St = gmshCircleWithHole(C2,{C1,P2},0.05,[0.02,0.005],'gmsh_circle_with_circular_hole_point_tri');
Sq = gmshCircleWithHole(C2,{C1,P2},0.05,[0.02,0.005],'gmsh_circle_with_circular_hole_point_quad',2,'recombine');

figure('Name','Circle with circular hole and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with ellipsoidal hole and point
fprintf('\n2D Circle with ellipsoidal hole and point\n');
fprintf('\n');
St = gmshCircleWithHole(C2,{E1,P2},0.05,[0.02,0.005],'gmsh_circle_with_ellipsoidal_hole_point_tri');
Sq = gmshCircleWithHole(C2,{E1,P2},0.05,[0.02,0.005],'gmsh_circle_with_ellipsoidal_hole_point_quad',2,'recombine');

figure('Name','Circle with ellipsoidal hole and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with circular hole, crack and point
fprintf('\n2D Circle with circular hole, crack and point\n');
fprintf('\n');
St = gmshCircleWithHole(C2,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_circle_with_circular_hole_crack_point_tri');
Sq = gmshCircleWithHole(C2,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_circle_with_circular_hole_crack_point_quad',2,'recombine');

figure('Name','Circle with circular hole, crack and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with ellipsoidal hole, crack and point
fprintf('\n2D Circle with ellipsoidal hole, crack and point\n');
fprintf('\n');
St = gmshCircleWithHole(C2,{E1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_circle_with_ellipsoidal_hole_crack_point_tri');
Sq = gmshCircleWithHole(C2,{E1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_circle_with_ellipsoidal_hole_crack_point_quad',2,'recombine');

figure('Name','Circle with ellipsoidal hole, crack and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with circular hole, line and point
fprintf('\n2D Circle with circular hole, line and point\n');
fprintf('\n');
St = gmshCircleWithHole(C2,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_circle_with_circular_hole_line_point_tri',2,'noduplicate');
Sq = gmshCircleWithHole(C2,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_circle_with_circular_hole_line_point_quad',2,'noduplicate','recombine');

figure('Name','Circle with circular hole, line and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with ellipsoidal hole, line and point
fprintf('\n2D Circle with ellipsoidal hole, line and point\n');
fprintf('\n');
St = gmshCircleWithHole(C2,{E1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_circle_with_ellipsoidal_hole_crack_point_tri',2,'noduplicate');
Sq = gmshCircleWithHole(C2,{E1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_circle_with_ellipsoidal_hole_crack_point_quad',2,'noduplicate','recombine');

figure('Name','Circle with ellipsoidal hole, line and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with circular inclusion
fprintf('\n2D Circle with circular inclusion\n');
fprintf('\n');
St = gmshCircleWithInclusion(C2,C1,0.05,0.02,'gmsh_circle_with_circular_inclusion_tri');
Sq = gmshCircleWithInclusion(C2,C1,0.05,0.02,'gmsh_circle_with_circular_inclusion_quad',2,'recombine');

figure('Name','Circle with circular inclusion')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with ellipsoidal inclusion
fprintf('\n2D Circle with ellipsoidal inclusion\n');
fprintf('\n');
St = gmshCircleWithInclusion(C2,E1,0.05,0.02,'gmsh_circle_with_ellipsoidal_inclusion_tri');
Sq = gmshCircleWithInclusion(C2,E1,0.05,0.02,'gmsh_circle_with_ellipsoidal_inclusion_quad',2,'recombine');

figure('Name','Circle with ellipsoidal inclusion')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with circular inclusion and point
fprintf('\n2D Circle with circular inclusion and point\n');
fprintf('\n');
St = gmshCircleWithInclusion(C2,{C1,P2},0.05,[0.02,0.005],'gmsh_circle_with_circular_inclusion_point_tri');
Sq = gmshCircleWithInclusion(C2,{C1,P2},0.05,[0.02,0.005],'gmsh_circle_with_circular_inclusion_point_quad',2,'recombine');

figure('Name','Circle with circular inclusion and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with ellipsoidal inclusion and point
fprintf('\n2D Circle with ellipsoidal inclusion and point\n');
fprintf('\n');
St = gmshCircleWithInclusion(C2,{E1,P2},0.05,[0.02,0.005],'gmsh_circle_with_ellipsoidal_inclusion_point_tri');
Sq = gmshCircleWithInclusion(C2,{E1,P2},0.05,[0.02,0.005],'gmsh_circle_with_ellipsoidal_inclusion_point_quad',2,'recombine');

figure('Name','Circle with ellipsoidal inclusion and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with circular inclusion, line and point
fprintf('\n2D Circle with circular inclusion, line and point\n');
fprintf('\n');
St = gmshCircleWithInclusion(C2,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_circle_with_circular_inclusion_line_point_tri',2);
Sq = gmshCircleWithInclusion(C2,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_circle_with_circular_inclusion_line_point_quad',2,'recombine');

figure('Name','Circle with circular inclusion, line and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with ellipsoidal inclusion, line and point
fprintf('\n2D Circle with ellipsoidal inclusion, line and point\n');
fprintf('\n');
St = gmshCircleWithInclusion(C2,{E1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_circle_with_ellipsoidal_inclusion_line_point_tri');
Sq = gmshCircleWithInclusion(C2,{E1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_circle_with_ellipsoidal_inclusion_line_point_quad',2,'recombine');

figure('Name','Circle with ellipsoidal inclusion, line and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with single edge crack
fprintf('\n2D Domain with single edge crack\n');
fprintf('\n');
St = gmshDomainWithSingleEdgeCrack(D,L1,0.05,0.01,'gmsh_domain_with_single_edge_crack_tri');
Sq = gmshDomainWithSingleEdgeCrack(D,L1,0.05,0.01,'gmsh_domain_with_single_edge_crack_quad',2,'recombine');
Strefcrack = gmshDomainWithSingleEdgeCrack(D,L1,0.05,0.01,'gmsh_domain_with_single_edge_crack_refined_tri',2,'refinecrack');
Sqrefcrack = gmshDomainWithSingleEdgeCrack(D,L1,0.05,0.01,'gmsh_domain_with_single_edge_crack_refined_quad',2,'recombine','refinecrack');

figure('Name','Domain with single edge crack')
clf
subplot(2,2,1)
plotparamelem(St,'group')
subplot(2,2,2)
plotparamelem(Sq,'group')
subplot(2,2,3)
plotparamelem(Strefcrack,'group')
subplot(2,2,4)
plotparamelem(Sqrefcrack,'group')

%% Domain with single edge circular notch
fprintf('\n2D Domain with single edge circular notch\n');
fprintf('\n');
c = max(getsize(D))/100; % notch width
St = gmshDomainWithSingleEdgeNotch(D,L1,c,0.05,0.01,'gmsh_domain_with_single_edge_circular_notch_tri',2,'c');
Sq = gmshDomainWithSingleEdgeNotch(D,L1,c,0.05,0.01,'gmsh_domain_with_single_edge_circular_notch_quad',2,'c','recombine');
Strefcrack = gmshDomainWithSingleEdgeNotch(D,L1,c,0.05,0.01,'gmsh_domain_with_single_edge_circular_notch_refined_tri',2,'c','refinecrack');
Sqrefcrack = gmshDomainWithSingleEdgeNotch(D,L1,c,0.05,0.01,'gmsh_domain_with_single_edge_circular_notch_refined_quad',2,'c','recombine','refinecrack');

figure('Name','Domain with single edge circular notch')
clf
subplot(2,2,1)
plotparamelem(St,'group')
subplot(2,2,2)
plotparamelem(Sq,'group')
subplot(2,2,3)
plotparamelem(Strefcrack,'group')
subplot(2,2,4)
plotparamelem(Sqrefcrack,'group')

%% Domain with single edge rectangular notch
fprintf('\n2D Domain with single edge rectangular notch\n');
fprintf('\n');
c = max(getsize(D))/100; % notch width
St = gmshDomainWithSingleEdgeNotch(D,L1,c,0.05,0.01,'gmsh_domain_with_single_edge_rectangular_notch_tri',2,'r');
Sq = gmshDomainWithSingleEdgeNotch(D,L1,c,0.05,0.01,'gmsh_domain_with_single_edge_rectangular_notch_quad',2,'r','recombine');
Strefcrack = gmshDomainWithSingleEdgeNotch(D,L1,c,0.05,0.01,'gmsh_domain_with_single_edge_rectangular_notch_refined_tri',2,'r','refinecrack');
Sqrefcrack = gmshDomainWithSingleEdgeNotch(D,L1,c,0.05,0.01,'gmsh_domain_with_single_edge_rectangular_notch_refined_quad',2,'r','recombine','refinecrack');

figure('Name','Domain with single edge rectangular notch')
clf
subplot(2,2,1)
plotparamelem(St,'group')
subplot(2,2,2)
plotparamelem(Sq,'group')
subplot(2,2,3)
plotparamelem(Strefcrack,'group')
subplot(2,2,4)
plotparamelem(Sqrefcrack,'group')

%% Domain with single edge V notch
fprintf('\n2D Domain with single edge V notch\n');
fprintf('\n');
c = max(getsize(D))/100; % notch width
St = gmshDomainWithSingleEdgeNotch(D,L1,c,0.05,0.01,'gmsh_domain_with_single_edge_V_notch_tri',2,'v');
Sq = gmshDomainWithSingleEdgeNotch(D,L1,c,0.05,0.01,'gmsh_domain_with_single_edge_V_notch_quad',2,'v','recombine');
Strefcrack = gmshDomainWithSingleEdgeNotch(D,L1,c,0.05,0.01,'gmsh_domain_with_single_edge_V_notch_refined_tri',2,'v','refinecrack');
Sqrefcrack = gmshDomainWithSingleEdgeNotch(D,L1,c,0.05,0.01,'gmsh_domain_with_single_edge_V_notch_refined_quad',2,'v','recombine','refinecrack');

figure('Name','Domain with single edge V notch')
clf
subplot(2,2,1)
plotparamelem(St,'group')
subplot(2,2,2)
plotparamelem(Sq,'group')
subplot(2,2,3)
plotparamelem(Strefcrack,'group')
subplot(2,2,4)
plotparamelem(Sqrefcrack,'group')

%% Domain with double edge crack (two symmetric or asymmetric edge cracks)
fprintf('\n2D Domain with double edge crack (two asymmetric edge cracks)\n');
fprintf('\n');
St = gmshDomainWithDoubleEdgeCrack(D,La,Lb,0.05,0.01,'gmsh_domain_with_double_edge_crack_tri');
Sq = gmshDomainWithDoubleEdgeCrack(D,La,Lb,0.05,0.01,'gmsh_domain_with_double_edge_crack_quad',2,'recombine');
Strefcrack = gmshDomainWithDoubleEdgeCrack(D,La,Lb,0.05,0.01,'gmsh_domain_with_double_edge_crack_refined_tri',2,'refinecrack');
Sqrefcrack = gmshDomainWithDoubleEdgeCrack(D,La,Lb,0.05,0.01,'gmsh_domain_with_double_edge_crack_refined_quad',2,'recombine','refinecrack');

figure('Name','Domain with double edge crack (two asymmetric edge cracks)')
clf
subplot(2,2,1)
plotparamelem(St,'group')
subplot(2,2,2)
plotparamelem(Sq,'group')
subplot(2,2,3)
plotparamelem(Strefcrack,'group')
subplot(2,2,4)
plotparamelem(Sqrefcrack,'group')

%% Domain with double edge circular notch (two symmetric or asymmetric edge circular notches)
fprintf('\n2D Domain with double edge circular notch (two asymmetric edge circular notches)\n');
fprintf('\n');
c = max(getsize(D))/100; % notch width
St = gmshDomainWithDoubleEdgeNotch(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_double_edge_circular_notch_tri',2,'c');
Sq = gmshDomainWithDoubleEdgeNotch(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_double_edge_circular_notch_quad',2,'c','recombine');
Strefcrack = gmshDomainWithDoubleEdgeNotch(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_double_edge_circular_notch_refined_tri',2,'c','refinecrack');
Sqrefcrack = gmshDomainWithDoubleEdgeNotch(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_double_edge_circular_notch_refined_quad',2,'c','recombine','refinecrack');

figure('Name','Domain with double edge circular notch (two asymmetric edge circular notches)')
clf
subplot(2,2,1)
plotparamelem(St,'group')
subplot(2,2,2)
plotparamelem(Sq,'group')
subplot(2,2,3)
plotparamelem(Strefcrack,'group')
subplot(2,2,4)
plotparamelem(Sqrefcrack,'group')

%% Domain with double edge rectangular notch (two symmetric or asymmetric edge rectangular notches)
fprintf('\n2D Domain with double edge rectangular notch (two asymmetric edge rectangular notches)\n');
fprintf('\n');
c = max(getsize(D))/100; % notch width
St = gmshDomainWithDoubleEdgeNotch(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_double_edge_rectangular_notch_tri',2,'r');
Sq = gmshDomainWithDoubleEdgeNotch(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_double_edge_rectangular_notch_quad',2,'r','recombine');
Strefcrack = gmshDomainWithDoubleEdgeNotch(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_double_edge_rectangular_notch_refined_tri',2,'r','refinecrack');
Sqrefcrack = gmshDomainWithDoubleEdgeNotch(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_double_edge_rectangular_notch_refined_quad',2,'r','recombine','refinecrack');

figure('Name','Domain with double edge rectangular notch (two asymmetric edge rectangular notches)')
clf
subplot(2,2,1)
plotparamelem(St,'group')
subplot(2,2,2)
plotparamelem(Sq,'group')
subplot(2,2,3)
plotparamelem(Strefcrack,'group')
subplot(2,2,4)
plotparamelem(Sqrefcrack,'group')

%% Domain with double edge V notch (two symmetric or asymmetric edge V notches)
fprintf('\n2D Domain with double edge V notch (two asymmetric edge V notches)\n');
fprintf('\n');
c = max(getsize(D))/100; % notch width
St = gmshDomainWithDoubleEdgeNotch(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_double_edge_V_notch_tri',2,'v');
Sq = gmshDomainWithDoubleEdgeNotch(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_double_edge_V_notch_v_quad',2,'v','recombine');
Strefcrack = gmshDomainWithDoubleEdgeNotch(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_double_edge_V_notch_refined_tri',2,'v','refinecrack');
Sqrefcrack = gmshDomainWithDoubleEdgeNotch(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_double_edge_V_notch_refined_quad',2,'v','recombine','refinecrack');

figure('Name','Domain with double edge V notch (two asymmetric edge V notches)')
clf
subplot(2,2,1)
plotparamelem(St,'group')
subplot(2,2,2)
plotparamelem(Sq,'group')
subplot(2,2,3)
plotparamelem(Strefcrack,'group')
subplot(2,2,4)
plotparamelem(Sqrefcrack,'group')

%% Domain with interior cracks
fprintf('\n2D Domain with interior cracks\n');
fprintf('\n');
St = gmshDomainWithInteriorCrack(D,{L2,L3},0.05,0.01,'gmsh_domain_with_interior_cracks_tri');
Sq = gmshDomainWithInteriorCrack(D,{L2,L3},0.05,0.01,'gmsh_domain_with_interior_cracks_quad',2,'recombine');

figure('Name','Domain with interior cracks')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Asymmetric plate with single edge crack and three holes
fprintf('\n2D Asymmetric plate with single edge crack and three holes\n');
fprintf('\n');
unit = 1e-3; % [mm]
a = 1*unit; % crack length
b = 6*unit; % crack offset from the centerline
c = 0.025*unit; % notch width
clD = 0.2*unit; % characteristic length for domain
clC = c; % characteristic length for edge crack/notch
clH = c; % characteristic length for circular holes
St = gmshAsymmetricPlateWithSingleEdgeCrackThreeHoles(a,b,clD,clC,clH,unit,'gmsh_asymmetric_plate_with_single_edge_crack_three_holes_tri');
Sq = gmshAsymmetricPlateWithSingleEdgeCrackThreeHoles(a,b,clD,clC,clH,unit,'gmsh_asymmetric_plate_with_single_edge_crack_three_holes_quad',2,'recombine');

figure('Name','Asymmetric plate with single edge crack and three holes')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Asymmetric plate with single edge circular notch and three holes
fprintf('\n2D Asymmetric plate with single edge circular notch and three holes\n');
fprintf('\n');
St = gmshAsymmetricPlateWithSingleEdgeNotchThreeHoles(a,b,c,clD,clC,clH,unit,'gmsh_asymmetric_plate_with_single_edge_circular_notch_three_holes_tri',2,'c');
Sq = gmshAsymmetricPlateWithSingleEdgeNotchThreeHoles(a,b,c,clD,clC,clH,unit,'gmsh_asymmetric_plate_with_single_edge_circular_notch_three_holes_quad',2,'c','recombine');

figure('Name','Asymmetric plate with single edge circular notch and three holes')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Asymmetric plate with single edge rectangular notch and three holes
fprintf('\n2D Asymmetric plate with single edge rectangular notch and three holes\n');
fprintf('\n');
St = gmshAsymmetricPlateWithSingleEdgeNotchThreeHoles(a,b,c,clD,clC,clH,unit,'gmsh_asymmetric_plate_with_single_edge_rectangular_notch_three_holes_tri',2,'r');
Sq = gmshAsymmetricPlateWithSingleEdgeNotchThreeHoles(a,b,c,clD,clC,clH,unit,'gmsh_asymmetric_plate_with_single_edge_rectangular_notch_three_holes_quad',2,'r','recombine');

figure('Name','Asymmetric plate with single edge rectangular notch and three holes')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Asymmetric plate with single edge V notch and three holes
fprintf('\n2D Asymmetric plate with single edge V notch and three holes\n');
fprintf('\n');
St = gmshAsymmetricPlateWithSingleEdgeNotchThreeHoles(a,b,c,clD,clC,clH,unit,'gmsh_asymmetric_plate_with_single_edge_V_notch_three_holes_tri',2,'v');
Sq = gmshAsymmetricPlateWithSingleEdgeNotchThreeHoles(a,b,c,clD,clC,clH,unit,'gmsh_asymmetric_plate_with_single_edge_V_notch_three_holes_quad',2,'v','recombine');

figure('Name','Asymmetric plate with single edge V notch and three holes')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% L-shaped panel
fprintf('\n2D L-shaped panel\n');
fprintf('\n');
St = gmshLshape(0.05,'gmsh_L_shape_tri');
Sq = gmshLshape(0.05,'gmsh_L_shape_quad',2,'recombine');

figure('Name','L-shaped panel')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% L-shaped panel
fprintf('\n2D L-shaped panel\n');
fprintf('\n');
a = 250e-3; % half-length
b = 30e-3; % distance of applied load from the right edge
e = 100e-3; % thickness
St = gmshLshapedPanel(a,b,e,20e-3,5e-3,'gmsh_L_shaped_panel_tri');
Sq = gmshLshapedPanel(a,b,e,20e-3,5e-3,'gmsh_L_shaped_panel_quad',2,'recombine');

figure('Name','L-shaped panel')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Canister
fprintf('\n2D Canister\n');
fprintf('\n');
St = gmshCanister(0.02,0.04,0.02,0.01,'gmsh_canister_tri');
Sq = gmshCanister(0.02,0.04,0.02,0.01,'gmsh_canister_quad',2,'recombine');

figure('Name','Canister')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Multiple domains
fprintf('\n2D Multiple domains\n');
fprintf('\n');
St = gmshMulti({D1,D2,D3},0.02,'gmsh_multi_tri');
Sq = gmshMulti({D1,D2,D3},0.02,'gmsh_multi_quad',2,'recombine');

figure('Name','Multiple domains')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Canister with multiple domains
fprintf('\n2D Canister with multiple domains\n');
fprintf('\n');
St = gmshCanisterMulti({D1,D2,D3},0.02,0.04,0.02,0.01,0.02,'gmsh_canister_multi_tri');
Sq = gmshCanisterMulti({D1,D2,D3},0.02,0.04,0.02,0.01,0.02,'gmsh_canister_multi_quad',2,'recombine');

figure('Name','Canister with multiple inclusions')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')
