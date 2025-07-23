function varargout = gmshfile_occ(E,cl,numbersurface,varargin)
% function G = gmshfile_occ(E,cl,numbersurface)
% E : ELLIPSE
% cl : characteristic length

if nargin<=2 || isempty(numbersurface), numbersurface = 1; end

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
G = createdisk(G,[E.cx,E.cy,E.cz],E.a,numbersurface,E.b);
G = setmeshsize(G,cl);

if ischarin('recombine',varargin)
    G = recombinesurface(G,numbersurface);
end

varargout{1} = G;
varargout{2} = numbersurface;
