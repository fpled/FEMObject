function u = runfilemmg2d(u,ext,options)
% function u = runfilemmg2d(u,ext,options)

if nargin==1 || isempty(ext)
    ext = '.msh';
end
command = ['mmg2d_O3 ' getfile(u,ext)];

if nargin<3 || isempty(options)
    if strcmpi(ext,'.mesh') % .mesh input file (and .mesh output file)
        options = '-3dMedit 2'; % to load a 2D .mesh file created with Gmsh (coordinates in 3D) and to produce a Gmsh 2D .mesh file (coordinates in 3D)
    end
else
    if verLessThan('matlab','9.1') % compatibility (<R2016b)
        opt3dMedit = ~isempty(strfind(options,'-3dMedit'));
    else
        opt3dMedit = contains(options,'-3dMedit');
    end
    if ~opt3dMedit
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            meshOutput = ~isempty(strfind(options,'.mesh'));
        else
            meshOutput = contains(options,'.mesh');
        end
        if strcmpi(ext,'.mesh') % .mesh input file (and .mesh output file)
            options = [options ' -3dMedit 2']; % to load a 2D .mesh file created with Gmsh (coordinates in 3D) and to produce a Gmsh 2D .mesh file (coordinates in 3D)
        elseif strcmpi(ext,'.msh') && meshOutput % .msh input file and .mesh output file
            options = [options ' -3dMedit 1']; % to force mmg2d to produce a 2D .mesh file readable by Gmsh (not anymore compatible with Medit)
            % options = [options ' -3dMedit 2']; % it also works
        end
    end
end

command = [command ' ' options];

pathname = getfemobjectoptions('mmgpath');
command = fullfile(pathname,command);
dos(command);
u.ismesh = 1;
