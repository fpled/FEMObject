function [Ds,tau,randh] = lsrandomsplit_elem2D(ls,tol,varargin)
% function [Ds,tau,randh] = lsrandomsplit_elem2D(ls,tol,varargin)

D=getclassin('MODEL',varargin);
PCM = PCMODEL(RANDVARS(ls),'order',1,'pcg');
X = RANDVARS(ls);
lspc = project(ls,PCM);
lspc= lseval(lspc,D);
lspcval = getvalue(lspc);
U = RANDVARS();
for i=1:getM(X)
    U{i}=RVUNIFORM(0,1);
end
   eoutglob = [];
   einglob  = [];
   ecutglob = [];
   
for p=1:getnbgroupelem(D)
    connec = getconnec(D.groupelem{p});
    numelem = getnumber(D.groupelem{p});
    fprintf('Group element %d ',p)
    for e=1:size(connec,1)
        fprintf('%d/%d ',e,size(connec,1))
        nume = numelem(e);
        cone = connec(e,:);
        lspce = lspcval(cone);
        keyboard
I_elem.order = 0;
I_elem.way = zeros(0,2);

xi=calc_sommets(calc_xi(I_elem));
I_elem.signsommet = find_Iin(xi,U,X,lspcval);

I_elemtot = cell(0,1);
I_elemtot = recursivecut(I_elemtot,I_elem,tol,U,X,lspcval);

figure(2)
clf
for i=1:length(I_elemtot)
   xi = calc_sommets(calc_xi(I_elemtot{i})); 
   switch I_elemtot{i}.state
       case 'out'
       col = 'g';
       case 'in'
       col='y';    
       case 'cut'
       col='m';    
   end
   %patch('Faces',[1,2,3,4],'Vertices',xi,'FaceColor',col) 

end
%keyboard

    end


end
   
function Iin = find_Iin(xi,U,X,lspcval)
x1 = transfer(U{1},X{1},xi(:,1));
x2 = transfer(U{2},X{2},xi(:,2));
Iin = sign(full(double(randomeval(lspcval,[x1 x2],X))));
return

function [I_elemtot] = recursivecut(I_elemtot,I_elem,tol,varargin)

if ((1/2)^(I_elem.order)<=tol) || comp_I(I_elem.signsommet) 
    I = I_elem.signsommet(:,1);
    if sum(I==-1 | I==0) == length(I)
    I_elem.state = 'in';    
    elseif sum(I==1 | I==0) == length(I)
    I_elem.state = 'out';    
    else
    I_elem.state = 'cut';            
    end
        
    I_elemtot = [I_elemtot , {I_elem}];
else
    I_subelem = cut_I_elem(I_elem,varargin{:});
for i=1:length(I_subelem)
    I_elemtot=recursivecut(I_elemtot,I_subelem{i},tol,varargin{:});
end
end

%keyboard
return

function rep = comp_I(Iin)
rep=1;
for i=2:size(Iin,2)
rep = rep & all(Iin(:,i)==Iin(:,1))    ;
end    


return

function [xi] = calc_xi(I_elem)

xi = sum((I_elem.way-1).*repmat((1/2).^[1:I_elem.order]',1,size(I_elem.way,2)),1);
xi = [xi;xi+repmat((1/2)^(I_elem.order),1,size(I_elem.way,2))];

return

function [xi] = calc_sommets(xiuni)

xi = [xiuni([1;2;2;1],1),xiuni([1;1;2;2],2)];


return


function I_subelem = cut_I_elem(I_elem,varargin)

I_subelem{1}.order=I_elem.order+1;
I_subelem{1}.way=[I_elem.way;[1,1]];
I_subelem{2}.order=I_elem.order+1;
I_subelem{2}.way=[I_elem.way;[2,1]];
I_subelem{3}.order=I_elem.order+1;
I_subelem{3}.way=[I_elem.way;[2,2]];
I_subelem{4}.order=I_elem.order+1;
I_subelem{4}.way=[I_elem.way;[1,2]];
m= size(I_elem.way,2);
xisub = zeros(0,m);
for i=1:length(I_subelem)
    xisub=[xisub;calc_sommets(calc_xi(I_subelem{i}))];
end
Iin = find_Iin(xisub,varargin{:});

for i=1:length(I_subelem)
       I_subelem{i}.signsommet = Iin(:,(i-1)*2^m+[1:2^m]);
end    
    
return
   