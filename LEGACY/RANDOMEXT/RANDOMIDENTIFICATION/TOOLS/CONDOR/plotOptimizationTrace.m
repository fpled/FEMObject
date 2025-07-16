function plotOptimizationTrace(mytrace,fig)
%
% CONDOR OPTIMIZER v1.11
%
% Copyright (c) 2004, Frank Vanden Berghen
%
% License & help @ http://iridia.ulb.ac.be/~fvandenb/optimization/CONDORdownload.html
%
eval('fig;','fig=figure;');
figure(fig);
clf;
hold on;
%   dbstop plotOptimizationTrace 6;
nEvaluations=size(mytrace,1);
[fumin,iopt]=min(mytrace(:,3));
l=plot3(mytrace(:,1),mytrace(:,2),zeros(1,nEvaluations)+fumin);
set(l,'Color','k');
set(l,'LineStyle',':');
l=plot3(mytrace(:,1),mytrace(:,2),mytrace(:,3));
set(l,'LineWidth',2);
[xs,ys,zs]=sphere(10);
r=(max(mytrace(:,1))-min(mytrace(:,1)))*.02;
xs=xs*r;
r=(max(mytrace(:,2))-min(mytrace(:,2)))*.02;
ys=ys*r;
r=(max(mytrace(:,3))-fumin)*.025;
zs=zs*r;
l=surf(xs+mytrace(1,1),ys+mytrace(1,2),zs+mytrace(1,3));
set(l,'FaceColor','r');
set(l,'EdgeColor','none');
l=surf(xs+mytrace(iopt,1),ys+mytrace(iopt,2),zs+fumin);
set(l,'FaceColor','g');
set(l,'EdgeColor','none');
for i=1:nEvaluations
    l=line([mytrace(i,1) mytrace(i,1)],[mytrace(i,2) mytrace(i,2)],[mytrace(i,3) fumin]);
    set(l,'Color','k');
end
view(-60,70);
l=title(sprintf('View of the optimization path followed by the Condor Optimizer\n(red dot=starting point; green dot=optimal point)'));
set(l,'VerticalAlignment','top');
hold off;
drawnow;
