function M = gmsh2femobject(indim,file,dim,varargin)
% function M = gmsh2femobject(indim,file,dim,option)
% indim : dimension de l'espace de travail
% file : GMSHFILE
% dim : recuperation des elements de dimension dim
%
% option :
% 'subfaces' :  si on veut que les facets, ridges, vertices aient leurs
% propres facets, ridges, vertices.
% 'nofaces' : si on ne veut pas de facets, ridges, peaks (pour sauver de l'espace memoire)

if nargin==2
    dim = indim;
end

if ~iswritten(file) && ~ismesh(file)
    file = writefile(file);
    % file = openfile(file,'.geo');
end

if ~ismesh(file)
    file = mesh(file,max(dim),varargin{:});
    % file = openfile(file,'.msh');
end

M = gmsh2femobject_model(indim,getfile(file,'.msh'),dim,varargin{:});
