function f = nodalload(S,D,compo,v,varargin)
% function f = nodalload(S,D,compo,v)
% effort nodal
% D zone d'application de l'effort :
%      objet geometrique (POINT, LIGNE, PLAN, ...)
%      on applique l'effort sur les points de S contenu dans D
%      si il y a plusieurs noeuds on divise la valeur de la force par le nombre
%      de noeuds
% compo : nom des composantes de la force
% v : valeurs des composantes vecteur colonne
%     peut avoir plusieurs colonnes correspondant � diff�rents chargements

nofree = ischarin('nofree',varargin);

s = size(v);
if isa(D,'double')
    numn = D;
else
    numn = ispointin(D,POINT(S.node)); % noeuds ou est applique l'effort
end
repddl = findddl(S,compo,numn,'dual');
f = sparse(getnbddl(S),s(2));
f(repddl,:) = repmat(v/length(numn),length(numn),1);

% f = FEVECTOR(f,S);
if ~nofree
    f = freevector(S,f);
end

if length(numn)>1
    fprintf('force nodale repartie sur %3d noeuds',length(numn))
    fprintf('\n')
end
