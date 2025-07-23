function u = subsref(u,s)
% function u = subsref(u,s)

switch s.type
    case '.'
        u = get(u,s.subs);
    case {'()'};
        switch length(s.subs)
            case 1
                u.MYDOUBLEND = u.MYDOUBLEND(:,:,s.subs{1});
                
                if ~isa(s.subs{1},'char') && (size(squeeze(s.subs{1}),2)>1 || ndims(squeeze(s.subs{1}))>2)
                    % length(find(size(s.subs{1})~=1)
                    u.MYDOUBLEND = reshape(u.MYDOUBLEND,[size2D(u.MYDOUBLEND),size(s.subs{1})]);
                end
            otherwise
                u.MYDOUBLEND = u.MYDOUBLEND(:,:,s.subs{:});
        end
    otherwise
        error('non defini')
        
end

