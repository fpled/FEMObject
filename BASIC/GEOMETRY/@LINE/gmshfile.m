function G = gmshfile(L,cl,numberpoints,numberline)
% function G = gmshfile(L,cl,numberpoints,numberline)
% L : LINE
% cl : characteristic length

if nargin<3 || isempty(numberpoints), numberpoints = 1:2; end
if nargin<4 || isempty(numberline), numberline = 1; end

G = GMSHFILE();
P = getvertices(L);
G = createpoints(G,P(1:2),cl,numberpoints);
G = createline(G,numberpoints,numberline);
