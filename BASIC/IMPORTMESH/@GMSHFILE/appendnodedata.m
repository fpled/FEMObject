function u = appendnodedata(u,q,metrictype,file)
% function u = appendnodedata(u,q,metrictype,file)

if nargin<3
    metrictype = 1; % 1-component (scalar) field: the size map
end

if nargin==4
    u = setfile(u,file);
end
file = getfile(u,'.msh');

[~,name,~] = fileparts(file);

q = full(double(q(:)));
nbnodes = size(q,1);

fid = fopen(file,'a');
fprintf(fid,'$NodeData\n');
fprintf(fid,'1\n');                % number-of-string-tags (1 string tag)
fprintf(fid,'"%s:metric"\n',name);   % string-tag: the name of the view (that must contains the :metric keyword) 
fprintf(fid,'1\n');                % number-of-real-tags (1 real tag)
fprintf(fid,'0.0\n');                % real-tag: the time value (0.0)
fprintf(fid,'3\n');                % number-of-integer-tags (3 integer tags)
fprintf(fid,'0\n');                  % integer-tag: the time step index (starting at 0)
fprintf(fid,'%u\n',metrictype);      % integer-tag: the number of field components of the data in the view
fprintf(fid,'%u\n',nbnodes);         % integer-tag: the number of entities (nodes) in the view
fprintf(fid,'%u %f\n',[(1:nbnodes)',q]'); % real numbers giving the values associated with each node (node-number value)
fprintf(fid,'$EndNodeData');
fclose(fid);
