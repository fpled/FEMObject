function varargout = gmshfile_occ(S,cl,numbervolume,varargin)
% function G = gmshfile_occ(S,cl,numbervolume)
% S : SPHERE
% cl : characteristic length

if nargin<=2 || isempty(numbervolume), numbervolume = 1; end

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
G = createsphere(G,[S.cx,S.cy,S.cz],S.r,numbervolume,varargin{:});
G = setmeshsize(G,cl);

if ischarin('recombine',varargin)
    G = recombinesurface(G);
end

varargout{1} = G;
varargout{2} = numbervolume;
