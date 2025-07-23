function varargout = gmshfile_occ(C,cl,numbersurface,varargin)
% function G = gmshfile_occ(C,cl,numbersurface)
% C : CIRCLE
% cl : characteristic length

if nargin<=2 || isempty(numbersurface), numbersurface = 1; end

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
G = createdisk(G,[C.cx,C.cy,C.cz],C.r,numbersurface);
G = setmeshsize(G,cl);

if ischarin('recombine',varargin)
    G = recombinesurface(G,numbersurface);
end

varargout{1} = G;
varargout{2} = numbersurface;
