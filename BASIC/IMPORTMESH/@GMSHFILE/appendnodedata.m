function u = appendnodedata(u,q,metrictype,file)
% function u = appendnodedata(u,q,metrictype,file)
% Append NodeData named "<meshname>:metric" in Gmsh file with MSH2 ASCII file format.

if nargin<3
    metrictype = 1; % type of metric = 1, the size map is defined as 1-component (scalar) field
end
if nargin==4
    u = setfile(u,file);
end
file = getfile(u,'.msh');

[~,name,~] = fileparts(file);
metricname = [name,':metric'];

q = full(double(q));
[n,m] = size(q);

fid = fopen(file,'a');
fprintf(fid,'$NodeData\n');       % For now, Mmg works only with metrics defined at mesh vertices
fprintf(fid,'1\n');               % number-of-string-tags (1 string tag for Mmg)
fprintf(fid,'"%s"\n',metricname);   % string-tag: name of the metric field (that must contains the :metric keyword suffix)
fprintf(fid,'1\n');               % number-of-real-tags (1 real tag for Mmg)
fprintf(fid,'0.0\n');               % real-tag: time value 0.0 (ignored value by Mmg)
fprintf(fid,'3\n');               % number-of-integer-tags (3 integer tags for Mmg)
fprintf(fid,'0\n');                 % integer-tag: time step index starting at 0 (ignored value by Mmg)
fprintf(fid,'%u\n',metrictype);     % integer-tag: type of data = number of field components (1:scalar, 3:vector, 9:tensor) of the data in the view
fprintf(fid,'%u\n',n);              % integer-tag: number of data (that must be the same number than the number of mesh nodes)
for i=1:n
    % Isotropic metric
    % fprintf(fid,'%u %e\n',i,q(i)); % real numbers giving the metric value associated with each node (node-number value)
    % Anisotropic flattened metric
    row = q(i,:);                 % real numbers giving the metric values associated with each node (node-number values)
    fprintf(fid,'%u %e',i,row(1));
    for j=2:m
        fprintf(fid,' %e',row(j));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'$EndNodeData\n');
fclose(fid);
