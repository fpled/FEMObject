function s = CARTESIAN3D(varargin)
% function s = CARTESIAN3D(varargin)

if nargin==0
    s.name = 'Coordonnees cartesiennes 3D';
    sp = SYSCOORD(3);
    s = class(s,'CARTESIAN3D',sp);
elseif nargin==1 && isa(varargin{1},'CARTESIAN3D')
    s = varargin{1};
else    
    s.name = 'Coordonnees cartesiennes 3D';
    if isa(varargin{1},'VECTEUR') && isa(varargin{2},'VECTEUR') && isa(varargin{3},'VECTEUR')
        base = [getcompo(varargin{1}),getcompo(varargin{2}),getcompo(varargin{3})];
    end
    axis = {'X','Y','Z'};
    sp = SYSCOORD(base,axis);
    s = class(s,'CARTESIAN3D',sp);
end
