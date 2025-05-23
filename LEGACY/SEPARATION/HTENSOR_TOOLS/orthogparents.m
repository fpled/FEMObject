function x = orthogparents(x,ii)
% function x = orthogparents(x,t)

if(x.is_orthog == true)
    return;
end

x_parent = x.parent;
x_is_leaf = x.is_leaf;
x_is_left = x.is_left;


while ii ~= 1;
    % Calculate QR decomposition of U{ii} or matricized B{ii}, set
    % U{ii} or B{ii} to Q.
    
    if(x_is_leaf(ii))
        
        % Calculate QR-decomposition
        [Q, R] = qr(x.U{ii}, 0);
        
        % Make sure that rank does not become zero
        if(size(R, 1) == 0)
            Q = ones(size(Q, 1), 1);
            R = ones(1, size(R, 2));
        end
        
        % Set U{ii} to Q
        x.U{ii} = Q;
        
    else
        % Matricize B{ii}
        B_mat = matricize(x.B{ii}, [1 2], 3, false);
        
        % Compute QR decomposition
        [Q, R] = qr(B_mat, 0);
        
        % Adjust sizes
        tsize_new = size(x.B{ii});
        tsize_new(3) = size(Q, 2);
        
        % Reshape Q into tensor B{ii}
        x.B{ii} = dematricize(Q, tsize_new, [1 2], 3, false);
    end
    
    R = full(R);
    
    % Index of parent node
    ii_par = x_parent(ii);
    
    % Multiply R into the transfer tensor of the parent node.
    if(x_is_left(ii))
        % left child of parent node
        x.B{ii_par} = ttm(x.B{ii_par}, R, 1);
    else
        % right child of parent node
        x.B{ii_par} = ttm(x.B{ii_par}, R, 2);
    end
    ii = ii_par;
end

end