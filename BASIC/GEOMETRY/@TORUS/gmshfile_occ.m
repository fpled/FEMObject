function varargout = gmshfile_occ(T,cl,numbervolume,varargin)
% function G = gmshfile_occ(T,cl,numbervolume)
% T : TORUS
% cl : characteristic length

if nargin<=2 || isempty(numbervolume), numbervolume = 1; end

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
G = createtorus(G,[T.cx,T.cy,T.cz],[T.r1,T.r2],numbervolume,T.angle);
G = setmeshsize(G,cl);

if ischarin('recombine',varargin)
    G = recombinesurface(G);
end

varargout{1} = G;
varargout{2} = numbervolume;
