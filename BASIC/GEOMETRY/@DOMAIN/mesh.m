function [S,X,Y,Z] = mesh(B,varargin)
% function S = mesh(B,n1)
% function S = mesh(B,n1,n2)
% function S = mesh(B,n1,n2,n3)
% function S = mesh(B,n1,n2,n3,'indim',indim,...)
% ni : number of elements in dimension i
% indim : space dimension (B.indim by default)
% optional arguments for addelem : material, param, option

indim = getcharin('indim',varargin,B.indim);
varargin = delcharin('indim',varargin);

switch indim
    case 1
        S = MODEL('UNID');
    case 2
        S = MODEL('PLANE');
    case 3
        S = MODEL('TRID');
    otherwise
        error('Wrong space dimension')
end

switch B.dim
    case 1
        n1 = varargin{1};
        
        p1 = getcoord(B.P1);
        p2 = getcoord(B.P2);
        
        p = double([p1(1),p2(1)]);
        
        if length(n1)>1
            px = n1;
            n1 = length(px)-1;
            px = p(1) + (px-px(1))/(px(end)-px(1)) * (p(2)-p(1));
        else
            px = linspace(p(1),p(2),n1+1);
        end
        X = full(px(:)).';
        
        node = NODE(X(:),1:numel(X));
        
        elem = [1:n1;2:n1+1]';
        elem = int32(elem);
        
        S = addnode(S,node);
        S = addelem(S,'SEG2',elem,varargin{:});
        S = setparamgroupelem(S,'partition',0,1);
        
    case 2
        n1 = varargin{1};
        n2 = varargin{2};
        
        p1 = getcoord(B.P1);
        p2 = getcoord(B.P2);
        
        p = double([p1(1),p2(1),p1(2),p2(2)]);
        
        if length(n1)>1
            px = n1;
            n1 = length(px)-1;
            px = p(1) + (px-px(1))/(px(end)-px(1)) * (p(2)-p(1));
        else
            px = linspace(p(1),p(2),n1+1);
        end
        
        if length(n2)>1
            py = n2;
            n2 = length(py)-1;
            py = p(3) + (py-py(1))/(py(end)-py(1)) * (p(4)-p(3));
        else
            py = linspace(p(3),p(4),n2+1);
        end
        
        [X,Y] = meshgrid(px,py);
        
        node = NODE([X(:),Y(:)],1:numel(X));
        
        elem = zeros(n1*n2,4,'int32');
        n=0;
        for j=1:n1
            i=(1:n2)';
            connec = [(j-1)*(n2+1)+i,...
                (j)*(n2+1)+i,...
                (j)*(n2+1)+i+1,...
                (j-1)*(n2+1)+i+1]    ;
            elem(n+(1:n2),:) = connec;
            n=n+n2;
        end
        
        S = addnode(S,node);
        S = addelem(S,'QUA4',elem,varargin{:});
        S = setparamgroupelem(S,'partition',0,1);
        
    case 3
        n1 = varargin{1};
        n2 = varargin{2};
        n3 = varargin{3};
        
        p1 = getcoord(B.P1);
        p2 = getcoord(B.P2);
        
        p = double([p1(1),p2(1),p1(2),p2(2),p1(3),p2(3)]);
        
        if length(n1)>1
            px = n1;
            n1 = length(px)-1;
            px = p(1) + (px-px(1))/(px(end)-px(1)) * (p(2)-p(1));
        else
            px = linspace(p(1),p(2),n1+1);
        end
        
        if length(n2)>1
            py = n2;
            n2 = length(py)-1;
            py = p(3) + (py-py(1))/(py(end)-py(1)) * (p(4)-p(3));
        else
            py = linspace(p(3),p(4),n2+1);
        end
        
        if length(n3)>1
            pz = n3;
            n3 = length(pz)-1;
            pz = p(5) + (pz-pz(1))/(pz(end)-pz(1)) * (p(6)-p(5));
        else
            pz = linspace(p(5),p(6),n3+1);
        end
        
        [X,Y,Z] = meshgrid(px,py,pz);
        
        node = NODE([X(:),Y(:),Z(:)],1:numel(X));
        
        elem = zeros(0,8,'int32');
        for k=1:n3
            for j=1:n1
                i=(1:n2)';
                connec = [(k-1)*(n2+1)*(n1+1)+(j-1)*(n2+1)+i,...
                    (k-1)*(n2+1)*(n1+1)+(j)*(n2+1)+i,...
                    (k-1)*(n2+1)*(n1+1)+(j)*(n2+1)+i+1,...
                    (k-1)*(n2+1)*(n1+1)+(j-1)*(n2+1)+i+1,...
                    (k)*(n2+1)*(n1+1)+(j-1)*(n2+1)+i,...
                    (k)*(n2+1)*(n1+1)+(j)*(n2+1)+i,...
                    (k)*(n2+1)*(n1+1)+(j)*(n2+1)+i+1,...
                    (k)*(n2+1)*(n1+1)+(j-1)*(n2+1)+i+1] ;
                elem = [elem;connec];
                
            end
        end
        
        S = addnode(S,node);
        S = addelem(S,'CUB8',elem,varargin{:});
        S = setparamgroupelem(S,'partition',0,1);
end
