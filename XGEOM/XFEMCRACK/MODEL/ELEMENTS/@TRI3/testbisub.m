function testbisub(t,ls1,ls2)

[cinin,cinout,coutin,coutout,xlnodeplus]=bilsdivide_oneelem(ls1,ls2);%
nodecoord = [nodelocalcoord(t);xlnodeplus];

t = TRI3(NODE(nodelocalcoord(t)));
figure(1)
clf
subplot(3,2,1)
patch('Faces',cinin,'Vertices',nodecoord,'FaceColor','r','EdgeColor','k');
plot(t,NODE(nodelocalcoord(t)))
axis square
title('inin')
subplot(3,2,2)
patch('Faces',cinout,'Vertices',nodecoord,'FaceColor','y','EdgeColor','k');
plot(t,NODE(nodelocalcoord(t)))
axis square
title('inout')
subplot(3,2,3)
patch('Faces',coutin,'Vertices',nodecoord,'FaceColor','g','EdgeColor','k');
plot(t,NODE(nodelocalcoord(t)))
axis square
title('outin')

subplot(3,2,4)
patch('Faces',coutout,'Vertices',nodecoord,'FaceColor','m','EdgeColor','k');
plot(t,NODE(nodelocalcoord(t)))
axis square
title('outout')

[ginin,ginout,goutin,goutout]=calc_bilssubgauss(t,ls1,ls2,1);

subplot(3,2,5:6)
patch('Faces',cinin,'Vertices',nodecoord,'FaceColor','r','EdgeColor','k');
plot(POINT(ginin.coord),'yo');
patch('Faces',cinout,'Vertices',nodecoord,'FaceColor','y','EdgeColor','k');
plot(POINT(ginout.coord),'ro');
patch('Faces',coutin,'Vertices',nodecoord,'FaceColor','g','EdgeColor','k');
plot(POINT(goutin.coord),'bo');
patch('Faces',coutout,'Vertices',nodecoord,'FaceColor','m','EdgeColor','k');
plot(POINT(goutout.coord),'ko');

axis square
