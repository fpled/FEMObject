function [se,fieldtype,fieldstorage,fieldddl] = sigma(elem,node,q,varargin)
% function [se,fieldtype,fieldstorage,fieldddl] = sigma(elem,node,q,varargin)

fieldtype = 'ddlgaussdual';
fieldddl = getddlgaussdual(elem);
if ischarin('node',varargin)
    gauss.coord = permute(nodelocalcoord(elem),[4,2,3,1]);
    gauss.nbgauss = getnbnode(elem);
    fieldstorage = 'node';
elseif ischarin('smooth',varargin)
    n = getcharin('intorder',varargin,'mass');
    gauss = calc_gauss(elem,n);
    fieldstorage = 'node';
else
    n = getcharin('intorder',varargin,'rigi');
    gauss = calc_gauss(elem,n);
    fieldstorage = 'gauss';
end

xnode = getcoord(node,getconnec(elem)');

qe = localize(elem,q);
mat = getmaterial(elem);
se = sigma(mat,elem,xnode,gauss.coord,qe,varargin{:});

if strcmp(getlstype(elem),'out')
    se(:)=0;
end
if ischarin('smooth',varargin)
    se = smooth(elem,node,se,gauss,'type',fieldtype);
end
