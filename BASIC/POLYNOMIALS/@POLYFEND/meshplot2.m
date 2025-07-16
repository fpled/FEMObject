function meshplot2(H,varargin)

matdecoup = getcharin('matdecoup',varargin);
varargin = delonlycharin('matdecoup',varargin);
stodim = getM(H);
vectdecoup = setdiff(1:stodim,matdecoup(:,1));
if ischarin('state',varargin)
            varargin = delonlycharin('state',varargin);
            state=1;
        else 
            state=0;
end
hold on
nbelem = length(H.e);
stateH = getstate(H);
% maxorder = zeros(1,stodim);
% for i = 1:nbelem
%     maxorder = max(maxorder,H.e{i}.order(end,:));
% end
    
for i =1:nbelem
%    test = zeros(size(matdecoup(:,2))); 
   elem = H.e{i};
   maxorder = elem.order; 
   decoup = 0:1/2^(maxorder(matdecoup(1,1))):1;
   test = decoup(find(decoup>=matdecoup(1,2) & decoup<=matdecoup(1,2)+1/2^(maxorder(matdecoup(:,1)))));
   
   I = calc_multi_indices(2,1,2);
   I=I(:,1:end-1);
   xi = calc_xi_ani(elem,stodim);
%keyboard
   testxi2 = xi(2,matdecoup(:,1));
   if  all(test==testxi2)
       for j=1:2
               for k=1:4
                   if I(k,j)==0
                       I(k,j) = xi(1,vectdecoup(j));
                   elseif I(k,j)==1
                       I(k,j) = xi(2,vectdecoup(j));
                   end
               end
       end
       vertex = [I(1,:);I(2,:);I(4,:);I(3,:)];
       faces = [ 1 2 3 4 ; 2 3 4 1];
       if state==1
                %keyboard
                state_elem = stateH(i);
                switch state_elem
                    case -1
                    col = [1 1 0];
                    case 1
                    col = [0 1 0];
                    case 0
                    col = [1 0 1];
                end
                options = {'FaceColor' 'flat' 'FaceVertexCData' col};
                %options = [options , varargin];
            else
                options = {'FaceColor'  'none'};
                %options = [options , varargin];
       end
            patch('Faces',faces,'Vertices',vertex,options{:})
            keyboard
   end

end

function [xi] = calc_xi_ani(I_elem,stodim)

xi = sum((I_elem.way-1).*((1/2).^(I_elem.order)),1);
if size(I_elem.order,1)==0
xi = [xi;xi+(1/2).^zeros(1,stodim)];    
else
xi = [xi;xi+(1/2).^(I_elem.order(end,:))];
end

return

function [xi] = calc_sommets(xiuni,stodim)
I=calc_multi_indices(stodim,1,2);
I=I(:,1:end-1);
xi=[];
for i=1:stodim
xi = [xi,xiuni(I(:,i)+1,i)];
end

return