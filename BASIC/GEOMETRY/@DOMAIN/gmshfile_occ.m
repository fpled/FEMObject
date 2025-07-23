function varargout = gmshfile_occ(D,cl,number,varargin)
% function G = gmshfile_occ(D,cl,number)
% D : DOMAIN
% cl : characteristic length

if nargin<=2 || isempty(number), number = 1; end

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
P1 = getvertex(D,1);
if D.dim==2
    G = createrectangle(G,[P1,0],getsize(D),number,varargin{:});
elseif D.dim==3
    G = createbox(G,P1,getsize(D),number,varargin{:});
end
G = setmeshsize(G,cl);

if ischarin('recombine',varargin)
    if D.dim==2
        G = recombinesurface(G,number);
    elseif D.dim==3
        G = recombinesurface(G);
    end
end

varargout{1} = G;
varargout{2} = number;
