
S = cast2matlab_model('PLANE','E:\PROGRAMMES\CASTEM\mesh_inclu_rand_0.txt');
mat=ELAS_ISOT('E',1,'NU',.3);
S = setmaterial(S,mat);
S = final(S);

getnbddl(S)
figure(1)
clf
subplot(1,2,1)
plot(S)

P  =POINT([0,0;1,0;1,1;0,1]);

L1 = LINE(P(1),P(2));
L2 = LINE(P(2),P(3));
L3 = LINE(P(3),P(4));
L4 = LINE(P(4),P(1));

S = addcl(S,P(1),'U',0);
S = addcl(S,P(2),'UY',0);
typeload=3;
switch typeload
case 1
    f = surfload(S,L2,'FX',1);
    f = f + surfload(S,L4,'FX',-1);
case 2
    f = surfload(S,L1,'FY',-1);
    f = f + surfload(S,L3,'FY',1);
case 3
    f = surfload(S,L1,'FX',-1);
    f = f + surfload(S,L2,'FY',1);
    f = f + surfload(S,L3,'FX',1);
    f = f + surfload(S,L4,'FY',-1);
end
K = calc_rigi(S);

u = K\f;
s = calc_sigma(S,u);
figure(2)
%subplot(1,2,2)
colorbar
plot(s,S+u/20);




