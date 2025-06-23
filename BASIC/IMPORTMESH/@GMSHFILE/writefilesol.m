function u = writefilesol(u,indim,q,file)
% function u = writefilesol(u,indim,q,file)

if nargin==4
    u = setfile(u,file);
end
file = getfile(u,'.sol');

q = full(double(q));
[n,m] = size(q);

% meshfile = getfile(u,'.mesh');
% [dim,nVertices] = readMeditMeshInfo(meshfile);
% if indim~=dim
%     error('writefilesol: sol dimension=%d but mesh dimension=%d.',indim,dim);
% end
% if n~=nVertices
%     error('writefilesol: number of metrics=%d but number of mesh nodes=%d.',n,nVertices);
% end

% Infer solution type
if m==1
    soltype = 1;  % scalar
elseif m==indim*(indim+1)/2
    soltype = 3;  % symmetric tensor
elseif m==indim
    soltype = 2;  % vector (unlikely for metric, but we allow)
else
    error('writefilesol: Unexpected number of columns=%d in solution for indim=%d.',m,indim);
end

nprec = 2;

% Open and write .sol file
fid = fopen(file,'w');
if fid < 0
    error('Cannot open solution file "%s" for writing.',file);
end
fprintf(fid,'MeshVersionFormatted %u\n\n',nprec); % nprec=1: simple precision, nprec=2: double precision
fprintf(fid,'Dimension %u\n\n',indim);            % metric dimension: must match with the mesh dimension
fprintf(fid,'SolAtVertices\n');                   % For now, Mmg works only with metrics defined at mesh vertices
fprintf(fid,'%u\n',n);                            % number of metrics: must match with the number of nodes (vertices) in the mesh file
fprintf(fid,'1 %u\n',soltype);                    % number of solutions per node (always 1 for Mmg) and type of each solution (1:scalar, 2:vector, 3:tensor). Note that vectors can not be used to prescribe a size map.
for i=1:n
    % Isotropic metric
    % fprintf(fid,'%e\n',q(i));                     % metric value at each node
    % Anisotropic flattened metric
    row = q(i,:);                                 % metric values at each node
    fprintf(fid,'%e',row(1));
    for j=2:m
        fprintf(fid,' %e',row(j));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
fprintf(fid,'End\n');
fclose(fid);

% fprintf('Wrote solution file "%s": %d vertices, dimension=%d, soltype=%d.\n',file,n,indim,soltype);
end

function [dim,nVertices] = readMeditMeshInfo(meshfile)
% Read Dimension and vertex count from a Medit .mesh file.
if verLessThan('matlab','9.1') % compatibility (<R2016b)
    strstrip = @strtrim;
    isstart = @(str,pat) strncmpi(str,pat,numel(pat));
else
    strstrip = @strip;
    isstart = @(str,pat) startsWith(str,pat,'IgnoreCase',true);
end

fid = fopen(meshfile,'r');
if fid<0
    error('Cannot open mesh file "%s".',meshfile);
end
dim = [];
nVertices = [];
while true
    tline = fgetl(fid);
    if ~ischar(tline), break; end
    t = strstrip(tline);
    % Detect Dimension
    if isstart(t,'Dimension')
        % Try to parse number on same line:
        nums = sscanf(t,'Dimension %d');
        if ~isempty(nums)
            dim = nums(1);
        else
            % read next non-empty line for the number
            t2 = '';
            while true
                nxt = fgetl(fid);
                if ~ischar(nxt), break; end
                t2 = strstrip(nxt);
                if ~isempty(t2), break; end
            end
            if isempty(t2)
                fclose(fid);
                error('Unexpected EOF after "Dimension" in mesh file.');
            end
            nums2 = sscanf(t2, '%d');
            if isempty(nums2)
                fclose(fid);
                error('Invalid dimension value: "%s".',t2);
            end
            dim = nums2(1);
        end
        % Continue to look for Vertices unless both found
    end
     % Detect Vertices
     if isstart(t,'Vertices')
         % Next non-empty line is count
         t2 = '';
         while true
             nxt = fgetl(fid);
             if ~ischar(nxt), break; end
             t2 = strstrip(nxt);
             if ~isempty(t2), break; end
         end
         if isempty(t2)
             fclose(fid);
             error('Unexpected EOF after "Vertices" in mesh file.');
         end
         nv = sscanf(t2, '%d');
         if isempty(nv)
             fclose(fid);
             error('Invalid vertex count: "%s".',t2);
         end
         nVertices = nv(1);
         % If we already have dim, we can break; else continue until both found
         if ~isempty(dim)
             break;
         end
     end
end
fclose(fid);
if isempty(dim)
    error('Could not find "Dimension" in mesh file "%s".',meshfile);
end
if isempty(nVertices)
    error('Could not find "Vertices" section in mesh file "%s".',meshfile);
end
end

