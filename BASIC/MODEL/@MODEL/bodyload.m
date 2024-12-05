function f = bodyload(S,D,compo,fun,varargin)
% function f = bodyload(S,D,compo,fun,varargin)
% assemblage d'une force volumique
% D zone d'application de l'effort :
%       objet geometrique de meme dimension que S
% compo : nom des composantes de la force
% fun : pointeur sur une fonction donnant la valeur des composantes en
%       colonnes (on peut avoir N colonnes correspondant a N cas de
%       chargement)
%       -> appel de fun(x,varargin{:}) ou x sont le tableau de coordonnees de P points
%       ou on veut evaluer la fonction
%       size(x,1)=P
%       size(x,2)=dimension de l'espace de definition des points
%       la sortie est un 3D array de taille length(compo)*N*P

nofree = ischarin('nofree',varargin);

if isa(fun,'FENODEFIELD')
    fun = FENODEFIELD(unfreevector(S,getvalue(fun)));
else
    if isa(fun,'FEELEMFIELD')
        varargin{1} = fun;
        nargin = 1;
        if size(double(fun),4)==1
            % fun = inline('repmat(double(a),[1 1 1 size(x,4) size(x,5)])','x','a');
            fun = @(x,a) repmat(double(a),[1 1 1 size(x,4) size(x,5)]);
        else
            % fun = inline('repmat(double(a),[1 1 1 1 size(x,5)])','x','a');
            fun = @(x,a) repmat(double(a),[1 1 1 1 size(x,5)]);
        end
    elseif isa(fun,'double')
        varargin{1} = fun;
        nargin = 1;
        % fun = inline('repmat(a,[1 1 size(x,3) size(x,4) size(x,5)])','x','a');
        fun = @(x,a) repmat(a,[1 1 size(x,3) size(x,4) size(x,5)]);
    end
    fun = fcnchk(fun); % transformation de la fonction en inline (fonction matlab)
end

if ~isempty(D)
    S0 = intersect(S,D); % partie du maillage ou on applique l'effort
    S0 = createddlnode(S0,S.node); % on copie les ddl de S
    P = calc_P(S,S0);
else
    S0 = S;
    P = 1;
end

f = calc_vector(S0,@load,compo,fun,varargin{:});

f = P'*f;

if ~nofree
    f = freevector(S,f);
end