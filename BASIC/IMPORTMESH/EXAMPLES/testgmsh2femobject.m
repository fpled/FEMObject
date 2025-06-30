%% Dimension 2
D = DOMAIN(2);
P1 = POINT([0.3,0.5]);
P2 = POINT([0.1,0.5]);
C1 = CIRCLE(0.5,0.5,.3);
C2 = CIRCLE(0.5,0.5,.5);
L1 = LIGNE([0.0,0.5],P1);
L2 = LIGNE([0.2,0.3],[0.6,0.7]);
L3 = LIGNE([0.7,0.2],[0.9,0.1]);

b = 1; % domain length
a = 0.2*b; % crack length
e = 0.025*b; % eccentricity
La = LIGNE([0.0,b/2+e],[a,b/2+e]);
Lb = LIGNE([b-a,b/2-e],[b,b/2-e]);

D1 = DOMAIN(2,[0.9,0.45],[1.0,0.55]);
D2 = DOMAIN(2,[0.51,0.45],[0.61,0.55]);
D3 = DOMAIN(2,[0.1,0.45],[0.2,0.55]);

%% Domain
St = gmsh(D,0.05,'filename','gmsh_domain_tri');
Sq = gmsh(D,0.05,'filename','gmsh_domain_quad','recombine');

figure('Name','Domain')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with point
St = gmsh(D,P1,0.05,0.01,'filename','gmsh_domain_with_point_tri');
Sq = gmsh(D,P1,0.05,0.01,'filename','gmsh_domain_with_point_quad','recombine');

figure('Name','Domain with point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle
St = gmsh(C1,0.05,'filename','gmsh_circle_tri');
Sq = gmsh(C1,0.05,'filename','gmsh_circle_quad','recombine');

figure('Name','Circle')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with point
St = gmsh(C1,P1,0.05,0.005,'filename','gmsh_circle_with_point');
Sq = gmsh(C1,P1,0.05,0.005,'filename','gmsh_circle_with_point_quad','recombine');

figure('Name','Circle with point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with hole
St = gmshDomainWithHole(D,C1,0.05,0.02,'gmsh_domain_with_hole_tri');
Sq = gmshDomainWithHole(D,C1,0.05,0.02,'gmsh_domain_with_hole_quad',2,'recombine');

figure('Name','Domain with hole')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with hole and point
St = gmshDomainWithHole(D,{C1,P2},0.05,[0.02,0.005],'gmsh_domain_with_hole_point_tri');
Sq = gmshDomainWithHole(D,{C1,P2},0.05,[0.02,0.005],'gmsh_domain_with_hole_point_quad',2,'recombine');

figure('Name','Domain with hole and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with hole, crack and point
St = gmshDomainWithHole(D,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_hole_crack_point_tri');
Sq = gmshDomainWithHole(D,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_hole_crack_point_quad',2,'recombine');

figure('Name','Domain with hole, crack and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with hole, line and point
St = gmshDomainWithHole(D,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_hole_line_point_tri',2,'noduplicate');
Sq = gmshDomainWithHole(D,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_hole_line_point_quad',2,'noduplicate','recombine');

figure('Name','Domain with hole, line and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with inclusion
St = gmshDomainWithInclusion(D,C1,0.05,0.02,'gmsh_domain_with_inclusion_tri');
Sq = gmshDomainWithInclusion(D,C1,0.05,0.02,'gmsh_domain_with_inclusion_quad',2,'recombine');

figure('Name','Domain with inclusion')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with inclusion and point
St = gmshDomainWithInclusion(D,{C1,P2},0.05,[0.02,0.005],'gmsh_domain_with_inclusion_point_tri');
Sq = gmshDomainWithInclusion(D,{C1,P2},0.05,[0.02,0.005],'gmsh_domain_with_inclusion_point_quad',2,'recombine');

figure('Name','Domain with inclusion and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with inclusion, line and point
St = gmshDomainWithInclusion(D,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_inclusion_line_point_tri');
Sq = gmshDomainWithInclusion(D,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_inclusion_line_point_quad',2,'recombine');

figure('Name','Domain with inclusion, line and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with hole
St = gmshCircleWithHole(C2,C1,0.05,0.02,'gmsh_circle_with_hole_tri');
Sq = gmshCircleWithHole(C2,C1,0.05,0.02,'gmsh_circle_with_hole_quad',2,'recombine');

figure('Name','Circle with hole')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with hole and point
St = gmshCircleWithHole(C2,{C1,P2},0.05,[0.02,0.005],'gmsh_circle_with_hole_point_tri');
Sq = gmshCircleWithHole(C2,{C1,P2},0.05,[0.02,0.005],'gmsh_circle_with_hole_point_quad',2,'recombine');

figure('Name','Circle with hole and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with inclusion
St = gmshCircleWithInclusion(C2,C1,0.05,0.02,'gmsh_circle_with_inclusion_tri');
Sq = gmshCircleWithInclusion(C2,C1,0.05,0.02,'gmsh_circle_with_inclusion_quad',2,'recombine');

figure('Name','Circle with inclusion')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with inclusion and point
St = gmshCircleWithInclusion(C2,{C1,P2},0.05,[0.02,0.005],'gmsh_circle_with_inclusion_point_tri');
Sq = gmshCircleWithInclusion(C2,{C1,P2},0.05,[0.02,0.005],'gmsh_circle_with_inclusion_point_quad',2,'recombine');

figure('Name','Circle with inclusion and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with single edge crack
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

%% Domain with two asymmetric edge cracks
St = gmshDomainWithTwoAsymmetricEdgeCracks(D,La,Lb,0.05,0.01,'gmsh_domain_with_two_asymmetric_edge_cracks_tri');
Sq = gmshDomainWithTwoAsymmetricEdgeCracks(D,La,Lb,0.05,0.01,'gmsh_domain_with_two_asymmetric_edge_cracks_quad',2,'recombine');
Strefcrack = gmshDomainWithTwoAsymmetricEdgeCracks(D,La,Lb,0.05,0.01,'gmsh_domain_with_two_asymmetric_edge_cracks_refined_tri',2,'refinecrack');
Sqrefcrack = gmshDomainWithTwoAsymmetricEdgeCracks(D,La,Lb,0.05,0.01,'gmsh_domain_with_two_asymmetric_edge_cracks_refined_quad',2,'recombine','refinecrack');

figure('Name','Domain with two asymmetric edge cracks')
clf
subplot(2,2,1)
plotparamelem(St,'group')
subplot(2,2,2)
plotparamelem(Sq,'group')
subplot(2,2,3)
plotparamelem(Strefcrack,'group')
subplot(2,2,4)
plotparamelem(Sqrefcrack,'group')

%% Domain with two asymmetric edge circular notches
c = max(getsize(D))/100; % notch width
St = gmshDomainWithTwoAsymmetricEdgeNotches(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_two_asymmetric_edge_circular_notches_tri',2,'c');
Sq = gmshDomainWithTwoAsymmetricEdgeNotches(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_two_asymmetric_edge_circular_notches_quad',2,'c','recombine');
Strefcrack = gmshDomainWithTwoAsymmetricEdgeNotches(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_two_asymmetric_edge_circular_notches_refined_tri',2,'c','refinecrack');
Sqrefcrack = gmshDomainWithTwoAsymmetricEdgeNotches(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_two_asymmetric_edge_circular_notches_refined_quad',2,'c','recombine','refinecrack');

figure('Name','Domain with two asymmetric edge circular notches')
clf
subplot(2,2,1)
plotparamelem(St,'group')
subplot(2,2,2)
plotparamelem(Sq,'group')
subplot(2,2,3)
plotparamelem(Strefcrack,'group')
subplot(2,2,4)
plotparamelem(Sqrefcrack,'group')

%% Domain with two asymmetric edge rectangular notches
c = max(getsize(D))/100; % notch width
St = gmshDomainWithTwoAsymmetricEdgeNotches(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_two_asymmetric_edge_rectangular_notches_tri',2,'r');
Sq = gmshDomainWithTwoAsymmetricEdgeNotches(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_two_asymmetric_edge_rectangular_notches_quad',2,'r','recombine');
Strefcrack = gmshDomainWithTwoAsymmetricEdgeNotches(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_two_asymmetric_edge_rectangular_notches_refined_tri',2,'r','refinecrack');
Sqrefcrack = gmshDomainWithTwoAsymmetricEdgeNotches(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_two_asymmetric_edge_rectangular_notches_refined_quad',2,'r','recombine','refinecrack');

figure('Name','Domain with two asymmetric edge rectangular notches')
clf
subplot(2,2,1)
plotparamelem(St,'group')
subplot(2,2,2)
plotparamelem(Sq,'group')
subplot(2,2,3)
plotparamelem(Strefcrack,'group')
subplot(2,2,4)
plotparamelem(Sqrefcrack,'group')

%% Domain with two asymmetric edge V notches
c = max(getsize(D))/100; % notch width
St = gmshDomainWithTwoAsymmetricEdgeNotches(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_two_asymmetric_edge_V_notches_tri',2,'v');
Sq = gmshDomainWithTwoAsymmetricEdgeNotches(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_two_asymmetric_edge_V_notches_v_quad',2,'v','recombine');
Strefcrack = gmshDomainWithTwoAsymmetricEdgeNotches(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_two_asymmetric_edge_V_notches_refined_tri',2,'v','refinecrack');
Sqrefcrack = gmshDomainWithTwoAsymmetricEdgeNotches(D,La,Lb,c,0.05,0.01,'gmsh_domain_with_two_asymmetric_edge_V_notches_refined_quad',2,'v','recombine','refinecrack');

figure('Name','Domain with two asymmetric edge V notches')
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
St = gmshDomainWithInteriorCrack(D,{L2,L3},0.05,0.01,'gmsh_domain_with_interior_cracks_tri');
Sq = gmshDomainWithInteriorCrack(D,{L2,L3},0.05,0.01,'gmsh_domain_with_interior_cracks_quad',2,'recombine');

figure('Name','Domain with interior cracks')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Asymmetric plate with single edge crack and three holes
unit = 1e-3; % [mm]
a = 1*unit; % crack length
b = 6*unit; % crack offset from the centerline
c = 0.025*unit; % notch width
clD = 0.2*unit; % characteristic length for domain
clC = c; % characteristic length for edge crack/notch
clH = c; % characteristic length for circular holes
St = gmshAsymmetricPlateWithSingleEdgeCrackThreeHoles(a,b,clD,clC,clH,unit,'gmsh_asymmetric_plate_with_single_edge_crack_three_holes_tri');
Sq = gmshAsymmetricPlateWithSingleEdgeCrackThreeHoles(a,b,clD,clC,clH,unit,'gmsh_asymmetric_plate_with_single_edge_crack_three_holes_quad',2,'recombine');

figure('Name','Asymmetric plate with single edge crack')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Asymmetric plate with single edge circular notch and three holes
St = gmshAsymmetricPlateWithSingleEdgeNotchThreeHoles(a,b,c,clD,clC,clH,unit,'gmsh_asymmetric_plate_with_single_edge_circular_notch_three_holes_tri',2,'c');
Sq = gmshAsymmetricPlateWithSingleEdgeNotchThreeHoles(a,b,c,clD,clC,clH,unit,'gmsh_asymmetric_plate_with_single_edge_circular_notch_three_holes_quad',2,'c','recombine');

figure('Name','Asymmetric plate with single edge circular notch and three holes')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Asymmetric plate with single edge rectangular notch and three holes
St = gmshAsymmetricPlateWithSingleEdgeNotchThreeHoles(a,b,c,clD,clC,clH,unit,'gmsh_asymmetric_plate_with_single_edge_rectangular_notch_three_holes_tri',2,'r');
Sq = gmshAsymmetricPlateWithSingleEdgeNotchThreeHoles(a,b,c,clD,clC,clH,unit,'gmsh_asymmetric_plate_with_single_edge_rectangular_notch_three_holes_quad',2,'r','recombine');

figure('Name','Asymmetric plate with single edge rectangular notch and three holes')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Asymmetric plate with single edge V notch and three holes
St = gmshAsymmetricPlateWithSingleEdgeNotchThreeHoles(a,b,c,clD,clC,clH,unit,'gmsh_asymmetric_plate_with_single_edge_V_notch_three_holes_tri',2,'v');
Sq = gmshAsymmetricPlateWithSingleEdgeNotchThreeHoles(a,b,c,clD,clC,clH,unit,'gmsh_asymmetric_plate_with_single_edge_V_notch_three_holes_quad',2,'v','recombine');

figure('Name','Asymmetric plate with single edge V notch and three holes')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% L-shaped panel
St = gmshLshape(0.05,'gmsh_L_shape_tri');
Sq = gmshLshape(0.05,'gmsh_L_shape_quad',2,'recombine');

figure('Name','L-shaped panel')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% L-shaped panel
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
St = gmshCanister(0.02,0.04,0.02,0.01,'gmsh_canister_tri');
Sq = gmshCanister(0.02,0.04,0.02,0.01,'gmsh_canister_quad',2,'recombine');

figure('Name','Canister')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Multiple domains
St = gmshMulti({D1,D2,D3},0.02,'gmsh_multi_tri');
Sq = gmshMulti({D1,D2,D3},0.02,'gmsh_multi_quad',2,'recombine');

figure('Name','Multiple domains')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Canister with multiple domains
St = gmshCanisterMulti({D1,D2,D3},0.02,0.04,0.02,0.01,0.02,'gmsh_canister_multi_tri');
Sq = gmshCanisterMulti({D1,D2,D3},0.02,0.04,0.02,0.01,0.02,'gmsh_canister_multi_quad',2,'recombine');

figure('Name','Canister with multiple inclusions')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Dimension 3
D = DOMAIN(3);
P = POINT([0.3,0.5,0.75]);
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
St = gmsh(D,0.1,'filename','gmsh_domain_tet');
% Sq = gmsh(D,0.1,'filename','gmsh_domain_cub','recombine');

figure('Name','Domain')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Domain with point
St = gmsh(D,P,0.1,0.01,'filename','gmsh_domain_with_point_tet');
% Sq = gmsh(D,P,0.1,0.01,'filename','gmsh_domain_with_point_cub','recombine');

figure('Name','Domain with point')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Domain with single edge crack
St = gmshDomainWithSingleEdgeCrack(D,Q,0.1,0.01,'gmsh_domain_with_single_edge_crack_tet');
% Sq = gmshDomainWithSingleEdgeCrack(D,Q,0.1,0.01,'gmsh_domain_with_single_edge_crack_cub',3,'recombine');
Strefcrack = gmshDomainWithSingleEdgeCrack(D,Q,0.1,0.01,'gmsh_domain_with_single_edge_crack_refined_tet',3,'refinecrack');
% Sqrefcrack = gmshDomainWithSingleEdgeCrack(D,Q,0.1,0.01,'gmsh_domain_with_single_edge_crack_refined_cub',3,'recombine','refinecrack');

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
c = max(getsize(D))/100; % notch width
St = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.01,'gmsh_domain_with_single_edge_circular_notch_tet',3,'c');
% Sq = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.01,'gmsh_domain_with_single_edge_circular_notch_cub',3,'c','recombine');
Strefcrack = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.01,'gmsh_domain_with_single_edge_circular_notch_refined_tet',3,'c','refinecrack');
% Sqrefcrack = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.01,'gmsh_domain_with_single_edge_circular_notch_refined_cub',3,'c','recombine','refinecrack');

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
c = max(getsize(D))/100; % notch width
St = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.01,'gmsh_domain_with_single_edge_rectangular_notch_tet',3,'r');
% Sq = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.01,'gmsh_domain_with_single_edge_rectangular_notch_cub',3,'r','recombine');
Strefcrack = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.01,'gmsh_domain_with_single_edge_rectangular_notch_refined_tet',3,'r','refinecrack');
% Sqrefcrack = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.01,'gmsh_domain_with_single_edge_rectangular_notch_refined_cub',3,'r','recombine','refinecrack');

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
c = max(getsize(D))/100; % notch width
St = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.01,'gmsh_domain_with_single_edge_V_notch_tet',3,'v');
% Sq = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.01,'gmsh_domain_with_single_edge_V_notch_cub',3,'v','recombine');
Strefcrack = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.01,'gmsh_domain_with_single_edge_V_notch_refined_tet',3,'v','refinecrack');
% Sqrefcrack = gmshDomainWithSingleEdgeNotch(D,Q,c,0.1,0.01,'gmsh_domain_with_single_edge_V_notch_refined_cub',3,'v','recombine','refinecrack');

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

%% Domain with two asymmetric edge cracks
St = gmshdomainwithtwoasymmetricedgecracks(D,Qa,Qb,0.1,0.01,'gmsh_domain_with_two_asymmetric_edge_cracks_tet');
% Sq = gmshdomainwithtwoasymmetricedgecracks(D,Qa,Qb,0.1,0.01,'gmsh_domain_with_two_asymmetric_edge_cracks_cub',3,'recombine');
Strefcrack = gmshdomainwithtwoasymmetricedgecracks(D,Qa,Qb,0.1,0.01,'gmsh_domain_with_two_asymmetric_edge_cracks_refined_tet',3,'refinecrack');
% Sqrefcrack = gmshdomainwithtwoasymmetricedgecracks(D,Qa,Qb,0.1,0.01,'gmsh_domain_with_two_asymmetric_edge_cracks_refined_cub',3,'recombine','refinecrack');

figure('Name','Domain with two asymmetric edge cracks')
clf
subplot(2,2,1)
plotparamelem(St,'group')
% subplot(2,2,2)
% plotparamelem(Sq,'group')
subplot(2,2,3)
plotparamelem(Strefcrack,'group')
% subplot(2,2,4)
% plotparamelem(Sqrefcrack,'group')

%% Domain with two asymmetric edge circular notches
c = max(getsize(D))/100; % notch width
St = gmshdomainwithtwoasymmetricedgenotches(D,Qa,Qb,c,0.1,0.01,'gmsh_domain_with_two_asymmetric_edge_circular_notches_tet',3,'c');
% Sq = gmshdomainwithtwoasymmetricedgenotches(D,Qa,Qb,c,0.1,0.01,'gmsh_domain_with_two_asymmetric_edge_circular_notches_cub',3,'c','recombine');
Strefcrack = gmshdomainwithtwoasymmetricedgenotches(D,Qa,Qb,c,0.1,0.01,'gmsh_domain_with_two_asymmetric_edge_circular_notches_refined_tet',3,'c','refinecrack');
% Sqrefcrack = gmshdomainwithtwoasymmetricedgenotches(D,Qa,Qb,c,0.1,0.01,'gmsh_domain_with_two_asymmetric_edge_circular_notches_refined_cub',3,'c','recombine','refinecrack');

figure('Name','Domain with two asymmetric edge circular notches')
clf
subplot(2,2,1)
plotparamelem(St,'group')
% subplot(2,2,2)
% plotparamelem(Sq,'group')
subplot(2,2,3)
plotparamelem(Strefcrack,'group')
% subplot(2,2,4)
% plotparamelem(Sqrefcrack,'group')

%% Domain with two asymmetric edge rectangular notches
c = max(getsize(D))/100; % notch width
St = gmshdomainwithtwoasymmetricedgenotches(D,Qa,Qb,c,0.1,0.01,'gmsh_domain_with_two_asymmetric_edge_rectangular_notches_tet',3,'r');
% Sq = gmshdomainwithtwoasymmetricedgenotches(D,Qa,Qb,c,0.1,0.01,'gmsh_domain_with_two_asymmetric_edge_rectangular_notches_cub',3,'r','recombine');
Strefcrack = gmshdomainwithtwoasymmetricedgenotches(D,Qa,Qb,c,0.1,0.01,'gmsh_domain_with_two_asymmetric_edge_rectangular_notches_refined_tet',3,'r','refinecrack');
% Sqrefcrack = gmshdomainwithtwoasymmetricedgenotches(D,Qa,Qb,c,0.1,0.01,'gmsh_domain_with_two_asymmetric_edge_rectangular_notches_refined_cub',3,'r','recombine','refinecrack');

figure('Name','Domain with two asymmetric edge rectangular notches')
clf
subplot(2,2,1)
plotparamelem(St,'group')
% subplot(2,2,2)
% plotparamelem(Sq,'group')
subplot(2,2,3)
plotparamelem(Strefcrack,'group')
% subplot(2,2,4)
% plotparamelem(Sqrefcrack,'group')

%% Domain with two asymmetric edge V notches
c = max(getsize(D))/100; % notch width
St = gmshdomainwithtwoasymmetricedgenotches(D,Qa,Qb,c,0.1,0.01,'gmsh_domain_with_two_asymmetric_edge_V_notches_tet',3,'v');
% Sq = gmshdomainwithtwoasymmetricedgenotches(D,Qa,Qb,c,0.1,0.01,'gmsh_domain_with_two_asymmetric_edge_V_notches_cub',3,'v','recombine');
Strefcrack = gmshdomainwithtwoasymmetricedgenotches(D,Qa,Qb,c,0.1,0.01,'gmsh_domain_with_two_asymmetric_edge_V_notches_refined_tet',3,'v','refinecrack');
% Sqrefcrack = gmshdomainwithtwoasymmetricedgenotches(D,Qa,Qb,c,0.1,0.01,'gmsh_domain_with_two_asymmetric_edge_V_notches_refined_cub',3,'v','recombine','refinecrack');

figure('Name','Domain with two asymmetric edge V notches')
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