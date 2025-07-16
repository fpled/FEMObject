function u = recombinesurface(u,number)
% function u = recombinesurface(u,number)

u = stepcounter(u);
if nargin==1 || isempty(number)
    s = 'Recombine Surface{:} ;\n';
else
    s = ['Recombine Surface{' num2str(number) '} ;\n'];
end
u = addstring(u,s);
