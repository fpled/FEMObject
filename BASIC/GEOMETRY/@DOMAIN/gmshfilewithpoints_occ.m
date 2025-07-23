function varargout = gmshfilewithpoints_occ(D,P,clD,clP,numberembeddedpoints,number,varargin)
% function G = gmshfilewithpoints_occ(D,P,clD,clP,numberembeddedpoints,number)
% D : DOMAIN
% P : POINT
% clD, clP : characteristic length

if ~iscell(P), P = {P}; end
if nargin<=3 || isempty(clP), clP = clD; end
if isscalar(clP), clP = repmat(clP,1,length(P)); end

if nargin<=4 || isempty(numberembeddedpoints)
    if D.dim==2
        numberembeddedpoints = 4+(1:length(P));
    elseif D.dim==3
        numberembeddedpoints = 8+(1:length(P));
    end
end
if nargin<=5 || isempty(number), number = 1; end

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
P1 = getvertex(D,1);
if D.dim==2
    G = createrectangle(G,[P1,0],getsize(D),number,varargin{:});
elseif D.dim==3
    G = createbox(G,P1,getsize(D),number,varargin{:});
end
G = setmeshsize(G,clD);

G = createpoints(G,P,clP,numberembeddedpoints);
if D.dim==2
    G = embedpointsinsurface(G,numberembeddedpoints,number);
elseif D.dim==3
    G = embedpointsinvolume(G,numberembeddedpoints,number);
end

if ischarin('recombine',varargin)
    if D.dim==2
        G = recombinesurface(G,number);
    elseif D.dim==3
        G = recombinesurface(G);
    end
end

varargout{1} = G;
varargout{2} = number;
