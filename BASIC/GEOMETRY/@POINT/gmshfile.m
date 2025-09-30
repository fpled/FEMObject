function G = gmshfile(P,cl,numberpoint)
% function G = gmshfile(P,cl,numberpoint)
% P : POINT
% cl : characteristic length

if nargin<3 || isempty(numberpoint), numberpoint = 1; end

G = GMSHFILE();
G = createpoint(G,P,cl,numberpoint);
