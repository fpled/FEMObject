%% Dimension 2
D = DOMAIN(2);
P1 = POINT([0.3,0.5]);
P2 = POINT([0.1,0.5]);
C1 = CIRCLE(0.5,0.5,.3);
C2 = CIRCLE(0.5,0.5,.5);
L1 = LIGNE([0.0,0.5],P1);
L2 = LIGNE([0.2,0.3],[0.6,0.7]);
L3 = LIGNE([0.7,0.2],[0.9,0.1]);

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
St = gmshdomainwithhole(D,C1,0.05,0.02,'gmsh_domain_with_hole_tri');
Sq = gmshdomainwithhole(D,C1,0.05,0.02,'gmsh_domain_with_hole_quad',2,'recombine');

figure('Name','Domain with hole')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with hole and point
St = gmshdomainwithhole(D,{C1,P2},0.05,[0.02,0.005],'gmsh_domain_with_hole_point_tri');
Sq = gmshdomainwithhole(D,{C1,P2},0.05,[0.02,0.005],'gmsh_domain_with_hole_point_quad',2,'recombine');

figure('Name','Domain with hole and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with hole, crack and point
St = gmshdomainwithhole(D,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_hole_crack_point_tri');
Sq = gmshdomainwithhole(D,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_hole_crack_point_quad',2,'recombine');

figure('Name','Domain with hole, crack and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with hole, line and point
St = gmshdomainwithhole(D,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_hole_line_point_tri',2,'noduplicate');
Sq = gmshdomainwithhole(D,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_hole_line_point_quad',2,'noduplicate','recombine');

figure('Name','Domain with hole, line and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with inclusion
St = gmshdomainwithinclusion(D,C1,0.05,0.02,'gmsh_domain_with_inclusion_tri');
Sq = gmshdomainwithinclusion(D,C1,0.05,0.02,'gmsh_domain_with_inclusion_quad',2,'recombine');

figure('Name','Domain with inclusion')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with inclusion and point
St = gmshdomainwithinclusion(D,{C1,P2},0.05,[0.02,0.005],'gmsh_domain_with_inclusion_point_tri');
Sq = gmshdomainwithinclusion(D,{C1,P2},0.05,[0.02,0.005],'gmsh_domain_with_inclusion_point_quad',2,'recombine');

figure('Name','Domain with inclusion and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with inclusion, line and point
St = gmshdomainwithinclusion(D,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_inclusion_line_point_tri');
Sq = gmshdomainwithinclusion(D,{C1,L3,P2},0.05,[0.02,0.01,0.005],'gmsh_domain_with_inclusion_line_point_quad',2,'recombine');

figure('Name','Domain with inclusion, line and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with hole
St = gmshcirclewithhole(C2,C1,0.05,0.02,'gmsh_circle_with_hole_tri');
Sq = gmshcirclewithhole(C2,C1,0.05,0.02,'gmsh_circle_with_hole_quad',2,'recombine');

figure('Name','Circle with hole')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with hole and point
St = gmshcirclewithhole(C2,{C1,P2},0.05,[0.02,0.005],'gmsh_circle_with_hole_point_tri');
Sq = gmshcirclewithhole(C2,{C1,P2},0.05,[0.02,0.005],'gmsh_circle_with_hole_point_quad',2,'recombine');

figure('Name','Circle with hole and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with inclusion
St = gmshcirclewithinclusion(C2,C1,0.05,0.02,'gmsh_circle_with_inclusion_tri');
Sq = gmshcirclewithinclusion(C2,C1,0.05,0.02,'gmsh_circle_with_inclusion_quad',2,'recombine');

figure('Name','Circle with inclusion')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Circle with inclusion and point
St = gmshcirclewithinclusion(C2,{C1,P2},0.05,[0.02,0.005],'gmsh_circle_with_inclusion_point_tri');
Sq = gmshcirclewithinclusion(C2,{C1,P2},0.05,[0.02,0.005],'gmsh_circle_with_inclusion_point_quad',2,'recombine');

figure('Name','Circle with inclusion and point')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with single edge crack
St = gmshdomainwithedgecrack(D,L1,0.05,0.01,'gmsh_domain_with_edge_crack_tri');
Sq = gmshdomainwithedgecrack(D,L1,0.05,0.01,'gmsh_domain_with_edge_crack_quad',2,'recombine');

figure('Name','Domain with single edge crack')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with single edge circular notch
c = max(getsize(D))/100; % notch width
St = gmshdomainwithedgenotch(D,L1,c,0.05,0.01,'gmsh_domain_with_edge_notch_c_tri',2,'c');
Sq = gmshdomainwithedgenotch(D,L1,c,0.05,0.01,'gmsh_domain_with_edge_notch_c_quad',2,'c','recombine');

figure('Name','Domain with single edge circular notch')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with single edge rectangular notch
c = max(getsize(D))/100; % notch width
St = gmshdomainwithedgenotch(D,L1,c,0.05,0.01,'gmsh_domain_with_edge_notch_r_tri',2,'r');
Sq = gmshdomainwithedgenotch(D,L1,c,0.05,0.01,'gmsh_domain_with_edge_notch_r_quad',2,'r','recombine');

figure('Name','Domain with single edge rectangular notch')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with single edge V notch
c = max(getsize(D))/100; % notch width
St = gmshdomainwithedgenotch(D,L1,c,0.05,0.01,'gmsh_domain_with_edge_notch_v_tri',2,'v');
Sq = gmshdomainwithedgenotch(D,L1,c,0.05,0.01,'gmsh_domain_with_edge_notch_v_quad',2,'v','recombine');

figure('Name','Domain with single edge V notch')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Domain with interior cracks
St = gmshdomainwithinteriorcrack(D,{L2,L3},0.05,0.01,'gmsh_domain_with_interior_cracks_tri');
Sq = gmshdomainwithinteriorcrack(D,{L2,L3},0.05,0.01,'gmsh_domain_with_interior_cracks_quad',2,'recombine');

figure('Name','Domain with interior cracks')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Asymmetric notched plate with single edge crack
unit = 1e-3; % [mm]
a = 1*unit; % crack length
b = 6*unit; % crack offset from the centerline
c = 0.025*unit; % notch width
clD = 0.2*unit; % characteristic length for domain
clC = c; % characteristic length for edge crack/notch
clH = c; % characteristic length for circular holes
St = gmshasymmetricnotchedplatewithedgecrack(a,b,clD,clC,clH,unit,'gmsh_asymmetric_notched_plate_with_edge_crack_tri');
Sq = gmshasymmetricnotchedplatewithedgecrack(a,b,clD,clC,clH,unit,'gmsh_asymmetric_notched_plate_with_edge_crack_quad',2,'recombine');

figure('Name','Asymmetric notched plate with single edge crack')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Asymmetric notched plate with single edge circular notch
St = gmshasymmetricnotchedplatewithedgenotch(a,b,c,clD,clC,clH,unit,'gmsh_asymmetric_notched_plate_with_edge_notch_c_tri',2,'c');
Sq = gmshasymmetricnotchedplatewithedgenotch(a,b,c,clD,clC,clH,unit,'gmsh_asymmetric_notched_plate_with_edge_notch_c_quad',2,'c','recombine');

figure('Name','Asymmetric notched plate with single edge circular notch')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Asymmetric notched plate with single edge rectangular notch
St = gmshasymmetricnotchedplatewithedgenotch(a,b,c,clD,clC,clH,unit,'gmsh_asymmetric_notched_plate_with_edge_notch_r_tri',2,'r');
Sq = gmshasymmetricnotchedplatewithedgenotch(a,b,c,clD,clC,clH,unit,'gmsh_asymmetric_notched_plate_with_edge_notch_r_quad',2,'r','recombine');

figure('Name','Asymmetric notched plate with single edge rectangular notch')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Asymmetric notched plate with single edge V notch
St = gmshasymmetricnotchedplatewithedgenotch(a,b,c,clD,clC,clH,unit,'gmsh_asymmetric_notched_plate_with_edge_notch_v_tri',2,'v');
Sq = gmshasymmetricnotchedplatewithedgenotch(a,b,c,clD,clC,clH,unit,'gmsh_asymmetric_notched_plate_with_edge_notch_v_quad',2,'v','recombine');

figure('Name','Asymmetric notched plate with single edge V notch')
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

%% Canister
St = gmshcanister(0.02,0.04,0.02,0.01,'gmsh_canister_tri');
Sq = gmshcanister(0.02,0.04,0.02,0.01,'gmsh_canister_quad',2,'recombine');

figure('Name','Canister')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Multiple domains
St = gmshmulti({D1,D2,D3},0.02,'gmsh_multi_tri');
Sq = gmshmulti({D1,D2,D3},0.02,'gmsh_multi_quad',2,'recombine');

figure('Name','Multiple domains')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Canister with multiple domains
St = gmshcanistermulti({D1,D2,D3},0.02,0.04,0.02,0.01,0.02,'gmsh_canister_multi_tri');
Sq = gmshcanistermulti({D1,D2,D3},0.02,0.04,0.02,0.01,0.02,'gmsh_canister_multi_quad',2,'recombine');

figure('Name','Canister with multiple inclusions')
clf
subplot(1,2,1)
plotparamelem(St,'group')
subplot(1,2,2)
plotparamelem(Sq,'group')

%% Dimension 3
D = DOMAIN(3);
P = POINT([0.3,0.5,0.75]);
a = 0.3;
b = 0.5;
e = 1;
Q = QUADRANGLE([0.0,b,0.0],[a,b,0.0],[a,b,e],[0.0,b,e]);
c = 0.01; % notch width

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
St = gmsh(D,P,0.1,0.05,'filename','gmsh_domain_with_point_tet');
% Sq = gmsh(D,P,0.1,0.05,'filename','gmsh_domain_with_point_cub','recombine');

figure('Name','Domain with point')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Domain with single edge crack
St = gmshdomainwithedgecrack(D,Q,0.1,0.005,'gmsh_domain_with_edge_crack_tet');
% Sq = gmshdomainwithedgecrack(D,Q,0.1,0.005,'gmsh_domain_with_edge_crack_cub',3,'recombine');

figure('Name','Domain with single edge crack')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Domain with single edge circular notch
St = gmshdomainwithedgenotch(D,Q,c,0.1,0.005,'gmsh_domain_with_edge_notch_c_tet',3,'c');
% Sq = gmshdomainwithedgenotch(D,Q,c,0.1,0.005,'gmsh_domain_with_edge_notch_c_cub',3,'c','recombine');

figure('Name','Domain with single edge circular notch')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Domain with single edge rectangular notch
St = gmshdomainwithedgenotch(D,Q,c,0.1,0.005,'gmsh_domain_with_edge_notch_r_tet',3,'r');
% Sq = gmshdomainwithedgenotch(D,Q,c,0.1,0.005,'gmsh_domain_with_edge_notch_r_cub',3,'r','recombine');

figure('Name','Domain with single edge rectangular notch')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')

%% Domain with single edge V notch
St = gmshdomainwithedgenotch(D,Q,c,0.1,0.005,'gmsh_domain_with_edge_notch_v_tet',3,'v');
% Sq = gmshdomainwithedgenotch(D,Q,c,0.1,0.005,'gmsh_domain_with_edge_notch_v_cub',3,'v','recombine');

figure('Name','Domain with single edge V notch')
clf
% subplot(1,2,1)
plotparamelem(St,'group')
% subplot(1,2,2)
% plotparamelem(Sq,'group')
