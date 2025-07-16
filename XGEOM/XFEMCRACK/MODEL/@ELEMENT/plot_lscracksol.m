function varargout = plot_lscracksol(elem,q,ls,varargin)
% function varargout = plot_lscracksol(elem,q,ls,varargin)

if ~isenrich(elem)
    Helem = plot_sol(elem,q,varargin{:})
else
    nodecoord=double(getcoord(node));
    connec = calc_conneclocal(elem,node) ;
    qe=localize(elem,q);
    xnode = node(elem);
    
    p = permute(nodelocalcoord(elem),[4,2,3,1]);
    u = double(calc_Nls(elem,xnode,p)*qe);
    u = reshape(u,[size(u,1),numel(u)/size(u,1)])';
    warning('mal programme')
    
    globconnec = reshape(1:numel(connec),size(connec));
    
    ampl = getcharin('ampl',varargin,1);
    varargin=delcharin('ampl',varargin);
    
    
    H = patch('Faces',globconnec,...
        'Vertices',nodecoord(connec(:),:)+ampl*u,varargin{:});
    
    if nargout ==1
        varargout{1}=H;
    end
    
end


