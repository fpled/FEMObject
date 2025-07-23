nbt = 10;
r=1;
Ld = 1;
Lt = 1/5;
bd = 0.1;
hd = 0.1;
Dt = 0.1;
mat=MATERIALS();
mat{1} = ELAS_BEAM('E',1,'NU',0.3,...
    'S',bd*hd,'IZ',bd*hd^3/12);
mat{2} = ELAS_ISOT('E',1,'NU',0.3,'S',pi*Dt^2/4);

P = POINT([0,0;Ld,0]);
rept = linspace(0,Ld,nbt+2)';
Pt = POINT([rept(2:end-1),ones(nbt,1)*(-Lt)]);

S1 = mesh(LINE(P(1),P(2)),(nbt+1)*r,mat{1});
S1 = convertelem(S1,'BEAM');
clear S2;
for k=1:nbt
    S2{k} = mesh(LINE(Pt(k),Pt(k)+VECTOR([0;Lt])),1,mat{2});
    S2{k} = convertelem(S2{k},'BARR'); 
end
S = union(S1,S2{:});

S = final(S);
S = addcl(S,Pt,'U');
S = addcl(S,P(1),'UX');

f = bodyload(S,S1,'FY',1);

K0 = calc_rigi(S,'selgroup',1);
K = K0 ; 

raid = RVLOGNORMAL(1,0.3,0,'stat');

RV = RANDVARS();
RV{1}=raid;


X = PCMODEL(RV,'order',4) ;  
PC = getPC(X);
% pour avoir une realisation de la matrice de rigidite associee a une
% realisation des rigidites (avec hypothese que toutes les 
% rigidites suivent la meme loi mais sont independantes)

% raidr=random(raid,nbt,1);
% clear Kt
% for i=1:nbt    
%    Kt{i} = calc_rigi(S,'selgroup',i+1);
%    K = K + Kt{i}*raidr(i);
% end

% pour avoir une realisation de la matrice de rigidite associee a une
% realisation des rigidites (avec hypothese que toutes les 
% rigidites sont une seule variable aleatoire)
% raidr=random(raid,1,1);
% Kt = calc_rigi(S,'selgroup',2:nbt+1);
% Kt = raidr*Kt;
% K = K0+Kt;

% resolution EF stochastique (avec hypothese que toutes les 
% rigidites sont une seule variable aleatoire)
Kt = calc_rigi(S,'selgroup',2:nbt+1);
K = K0 + X{1}*Kt;

% si on a les coeff sur le chaos (x=[1,4,...]) de la rigiditee (avec hypothese que toutes les 
% rigidites sont une seule variable aleatoire)
% k = PCMATRIX(x,[1,1],getPC(X));
% K = K0 + k*Kt;


if israndom(K)
    q = solve(K,f) ;
else
    q = K\f;
end

%% calcul des efforts dans les tirants
% ep = calc_epsilon(S,q);
% si = zeros(nbt,length(PC));
% for i=1:nbt
%     si{i} = double(ep{i+1}); % decomposition sur le chaos de la deformation dans le tirant i
% end

figure
clf
ampl=getsize(S)/max(abs(mean(q)));
plot(S,'Color','b')
plot(S+mean(q)*ampl,'Color','r')
