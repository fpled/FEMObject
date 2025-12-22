function S = build_model(D,varargin)
% function S = build_model(D,varargin)
% build MODEL S from GEOMOBJECT D
% function S = build_model(D,'nbelem',nbelem)
% function S = build_model(D,'nbelem',nbelem,'elemtype',elemtype,'indim',indim,varagin)
% function S = build_model(D,'cl',cl)
% function S = build_model(D,'cl',cl,'elemtype',elemtype,'indim',indim,'filename',filename,varargin)
% nbelem: number of elements
% cl: characteristic length
% elemtype: type of elements (optional, 'SEG2' in 1D, 'TRI3' in 2D, 'TET4' in 3D by default)
% filename: file name (optional, 'gmsh_file' by default)
% indim: space dimension (optional, getdim(D) by default)
% optional arguments for mesh : material, param, option

dim = getdim(D);
indim = getcharin('indim',varargin,getindim(D));
varargin = delcharin('indim',varargin);

switch dim
    case 0
        elemtype = getcharin('elemtype',varargin,'ELEMPOINT');
    case 1
        elemtype = getcharin('elemtype',varargin,'SEG2');
    case 2
        switch indim
            case 2
                elemtype = getcharin('elemtype',varargin,'TRI3');
            case 3
                elemtype = getcharin('elemtype',varargin,'DKT');
            otherwise
                error('Wrong space dimension')
        end
    case 3
        elemtype = getcharin('elemtype',varargin,'TET4');
    otherwise
        error('Wrong dimension')
end
varargin = delcharin('elemtype',varargin);

if ischarin('nbelem',varargin)
    nbelem = getcharin('nbelem',varargin);
    varargin = delcharin('nbelem',varargin);
    switch dim
        case 1
            S = mesh(D,nbelem(1),'indim',indim,varargin{:});
        case 2
            S = mesh(D,nbelem(1),nbelem(2),'indim',indim,varargin{:});
        case 3
            S = mesh(D,nbelem(1),nbelem(2),nbelem(3),'indim',indim,varargin{:});
    end
elseif ischarin('cl',varargin)
    cl = getcharin('cl',varargin);
    filename = getcharin('filename',varargin,'gmsh_file');
    varargin = delcharin({'cl','filename'},varargin);
    if ischarin('points',varargin)
        P = getcharin('points',varargin);
        varargin = delcharin('points',varargin);
        if ~strcmp(elemtype,'QUA4') && ~strcmp(elemtype,'CUB8') && ~strcmp(elemtype,'DKQ') && ~strcmp(elemtype,'DSQ') && ~strcmp(elemtype,'COQ4') && ~strcmp(elemtype,'STOKES')
            S = gmsh(D,P,cl,'filename',filename,'indim',indim,varargin{:});
        else
            S = gmsh(D,P,cl,'filename',filename,'indim',indim,'recombine',varargin{:});
        end
    else
        if ~strcmp(elemtype,'QUA4') && ~strcmp(elemtype,'CUB8') && ~strcmp(elemtype,'DKQ') && ~strcmp(elemtype,'DSQ') && ~strcmp(elemtype,'COQ4') && ~strcmp(elemtype,'STOKES')
            S = gmsh(D,cl,'filename',filename,'indim',indim,varargin{:});
        else
            S = gmsh(D,cl,'filename',filename,'indim',indim,'recombine',varargin{:});
        end
    end
    mat = getclassin('MATERIAL',varargin);
    option = getcharin('option',varargin,' ');
    S = setoption(S,option);
    S = setmaterial(S,mat);
else
    error('Wrong input arguments')
end

S = convertelem(S,elemtype,varargin{:});
