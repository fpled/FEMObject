function u = updatenodedata(u,q,metrictype,file)
% function u = updatenodedata(u,q,metrictype,file)
% Replace or append NodeData named "<meshname>:metric" in Gmsh file with MSH2 ASCII file format.

if verLessThan('matlab','9.1') % compatibility (<R2016b)
    strstrip = @strtrim;
    contain = @(str,pat) ~isempty(strfind(str,pat));
else
    strstrip = @strip;
    contain = @contains;
end

if nargin<3
    metrictype = 1; % type of metric = 1, the size map is defined as 1-component (scalar) field
end
if nargin==4
    u = setfile(u,file);
end
file = getfile(u,'.msh');

[~,name,~] = fileparts(file);
metricname = [name,':metric'];

q = full(double(q(:))); % ensure q is column vector double

% Open and read .msh file
fin = fopen(file,'r');
if fin<0
    error('Cannot open msh file "%s".',file);
end
tmp = [file,'.tmp'];
fout = fopen(tmp,'w');
if fout<0
    fclose(fin);
    error('Cannot open temp file.');
end

found = false;
while true
    tline = fgetl(fin);
    if ~ischar(tline)
        break;
    end

    if ~found && strcmp(strstrip(tline),'$NodeData')
        % Peek next two lines for the string tag
        line1 = fgetl(fin);
        line2 = fgetl(fin);
        if ~ischar(line1) || ~ischar(line2)
            % malformed: write back and continue
            fprintf(fout,'%s\n',tline);
            if ischar(line1), fprintf(fout,'%s\n',line1); end
            if ischar(line2), fprintf(fout,'%s\n',line2); end
            continue;
        end
        nstr = str2double(strstrip(line1));
        tagclean = strstrip(line2);
        if nstr==1
            % remove surrounding quotes if present
            if numel(tagclean)>=2 && ((tagclean(1)=='"' && tagclean(end)=='"') ...
                    || (tagclean(1)=='''' && tagclean(end)==''''))
                tagclean = tagclean(2:end-1);
            end
            if contain(tagclean,metricname)
                found = true;
                % skip rest of this block (1 real tag + real-tag line, 3 integer tags + integer-tag lines, data lines) until $EndNodeData
                % read real-tag count
                line3 = fgetl(fin);
                nreal = str2double(strstrip(line3));
                % skip real-tag line
                for k=1:nreal
                    fgetl(fin);
                end
                % read integer-tag count
                line4 = fgetl(fin);
                nint = str2double(strstrip(line4));
                % skip integer-tag lines
                ints = zeros(nint,1);
                for k=1:nint
                    il = fgetl(fin);
                    ints(k) = str2double(strstrip(il));
                end
                % read data count = last int
                if ~isempty(ints)
                    old_n = ints(end);
                else
                    old_n = 0;
                end
                % skip data lines
                for k=1:old_n
                    fgetl(fin);
                end
                % skip $EndNodeData
                fgetl(fin);
                % write new block
                writeBlock(fout, metricname, metrictype, q);
                continue;
            end
        end
        % If not matching, write header lines and copy until $EndNodeData
        fprintf(fout,'%s\n',tline);
        fprintf(fout,'%s\n',line1);
        fprintf(fout,'%s\n',line2);
        % copy until $EndNodeData
        while true
            l = fgetl(fin);
            if ~ischar(l), break; end
            fprintf(fout,'%s\n',l);
            if strcmp(strstrip(l),'$EndNodeData')
                break;
            end
        end
    else
        fprintf(fout,'%s\n',tline);
    end
end

% If not found, append at end
if ~found
    writeBlock(fout, metricname, metrictype, q);
end

fclose(fin);
fclose(fout);
movefile(tmp, file);
end

function writeBlock(fid,metricname,metrictype, q)
[n,m] = size(q);
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
end