function u = createcrack(u,dim,physicalgroup,openboundaryphysicalgroup,normalX,normalY,normalZ,newphysicalgroup,debugview,swaporientation,file)
% function u = createcrack(u,dim,physicalgroup,openboundaryphysicalgroup,normalX,normalY,normalZ,newphysicalgroup,debugview,swaporientation,file)

if nargin<2 || isempty(dim)
    dim = 1;
end
if nargin<3 || isempty(physicalgroup)
    physicalgroup = 1;
end
if nargin<4 || isempty(openboundaryphysicalgroup)
    openboundaryphysicalgroup = 0;
end
if nargin<5 || isempty(normalX)
    normalX = 0;
end
if nargin<6 || isempty(normalY)
    normalY = 0;
end
if nargin<7 || isempty(normalZ)
    normalZ = 1;
end
if nargin<8 || isempty(newphysicalgroup)
    newphysicalgroup = 0;
end
if nargin<9 || isempty(debugview)
    debugview = 0;
end
if nargin<10 || isempty(swaporientation)
    swaporientation = 0;
end
if nargin==11
    u = setfile(u,file);
end
% file = getfile(u,'.msh.opt');
file = getfile(u,'.geo','post');

fid = fopen(file,'w');
fprintf(fid,'Plugin(Crack).Dimension = %u;\n',dim);
fprintf(fid,'Plugin(Crack).PhysicalGroup = %u;\n',physicalgroup);
fprintf(fid,'Plugin(Crack).OpenBoundaryPhysicalGroup = %u;\n',openboundaryphysicalgroup);
fprintf(fid,'Plugin(Crack).NormalX = %f;\n',normalX);
fprintf(fid,'Plugin(Crack).NormalY = %f;\n',normalY);
fprintf(fid,'Plugin(Crack).NormalZ = %f;\n',normalZ);
fprintf(fid,'Plugin(Crack).NewPhysicalGroup = %u;\n',newphysicalgroup);
fprintf(fid,'Plugin(Crack).DebugView = %u;\n',debugview);
fprintf(fid,'Plugin(Crack).SwapOrientation = %u;\n',swaporientation);
fprintf(fid,'Plugin(Crack).Run;\n');
fclose(fid);

u.iswritten = 1;
