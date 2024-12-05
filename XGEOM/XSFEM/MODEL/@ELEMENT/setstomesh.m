function elem=setstomesh(elem,h,e)
% function elem=setstomesh(elem,h)
% elem : groupe d'�l�ments
% h : tableau de cellules , h{e} : POLYFEND correspondant au maillage
% stochastique de l'�l�ment e du groupe d'�l�ment 'elem'
if nargin==2
elem.stomesh=h;
elseif nargin ==3
    elem.stomesh{e}=h;
else
    error('mauvais arguments')
end