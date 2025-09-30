function varargout = gmsh(D,varargin)
% function u = gmsh(D,varargin)
% function u = gmsh(D,P,varargin)
% D : GEOMOBJECT
% P : POINT

filename    = getcharin('filename',varargin, 'gmsh_file');
indim       = getcharin('indim',varargin, getindim(D));
opencascade = any(ischarin({'occ','OpenCASCADE'},varargin));
extrude     = ischarin('extrude',varargin);
recombine   = ischarin('recombine',varargin);

varargin = delcharin({'filename','indim'},varargin);
varargin = delonlycharin({'extrude','recombine', ...
                          'occ','OpenCASCADE'},varargin);

flags = {};
if extrude,   flags{end+1} = 'extrude';   end
if recombine, flags{end+1} = 'recombine'; end

if ~isempty(varargin) && isscalar(varargin{1}) && ~isa(varargin{1},'POINT')
    % gmshfile
    if opencascade, fn = 'gmshfile_occ'; else, fn = 'gmshfile'; end
else
    % gmshfilewithpoints
    if opencascade, fn = 'gmshfilewithpoints_occ'; else, fn = 'gmshfilewithpoints'; end
end

if verLessThan('matlab','9.1') % compatibility (<R2016b)
    contain = @(str,pat) ~isempty(strfind(str,pat));
else
    contain = @contains;
end

% Locate the class‐specific .m file
allFiles = which([fn '.m'], '-all');
tag = ['@' class(D) '/' fn '.m'];
idx = find(contain(allFiles, tag), 1);
if isempty(idx)
    % fallback if no class‐folder match
    methodFile = fn;
else
    methodFile = allFiles{idx};
end

% Figure out how many positional inputs it expects
n_decl = nargin(methodFile);
if n_decl < 0
    n_pos = abs(n_decl) - 1; % has varargin
else
    n_pos = n_decl;          % fixed signature
end
n_pos = n_pos - 1; % subtract one for the object D itself
n_take = min(numel(varargin), n_pos);
fixed_args = varargin(1:n_take);
pads = repmat({[]}, 1, max(0, n_pos - n_take));

% Call the correct class‐method via feval
[G, number] = feval(fn, D, fixed_args{:}, pads{:}, flags{:});

switch D.dim
    case 2
        G = createphysicalsurface(G,number,1);
    case 3
        G = createphysicalvolume(G,number,1);
end

G = setfile(G, filename);
n = max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim, G, getdim(D):-1:getdim(D)-n+1);
 