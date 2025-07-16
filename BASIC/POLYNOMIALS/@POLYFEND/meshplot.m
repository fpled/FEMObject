function meshplot(H,varargin)
% meshplot(H,varargin)
% si la dimension stochastique est supérieure à 3 rentrer en varargin
% 'matdecoup' matrice de taille nsto-dimfixe*2 donnant les dimensions sto
% a fixer et les valeurs de ksi correspondantes pour la visualisation

plotpoint=ischarin('statepoint',varargin);
varargin = delonlycharin('statepoint',varargin);
state = ischarin('state',varargin);
varargin = delonlycharin('state',varargin);

if ~ischarin('matdecoup',varargin)
    switch getM(H)
        case 1
            I1 =getparam(getpoly(H,1),'I');
            j1 = H.j(:,1);
            I1 = I1(j1,:);
            node = I1;
            connec = zeros(H.n,2);
            connec(:) = [1 : numel(connec)];
            
        case 2
            I1 = getparam(getpoly(H,1),'I');
            I2 = getparam(getpoly(H,2),'I');
            j1 = H.j(:,1);
            j2 = H.j(:,2);
            I1 = I1(j1,:);
            I2 = I2(j2,:);
            node = [I1(:,1),I2(:,1);I1(:,2),I2(:,1);I1(:,2),I2(:,2);I1(:,1),I2(:,2)];
            connec = zeros(H.n,4);
            connec(:) = [1 : numel(connec)];
            
            if state
                col = getstate(H);
                col = col(:);
                options = {'FaceColor','flat','FaceVertexCData',col};
                options = [options , varargin];
            else
                options = {'FaceColor','none'};
                options = [options , varargin];
            end
            patch('Faces',connec,'Vertices',node,options{:})
            if plotpoint==1
                hold on
                Pcut2 = [];
                Pcut3 = [];
                Pcut4 = [];
                Pin = [];
                Pout = [];
                for k=1:H.n
                    xi1=calc_sommets(calc_xi_ani(H.e{k},2),2);
                    statepoint = H.statepoint(k,:);
                    for j=1:4
                        if statepoint(j)==2
                            Pcut2 = [Pcut2;xi1(j,:)];
                        elseif statepoint(j)==3
                            Pcut3 = [Pcut3;xi1(j,:)];
                        elseif statepoint(j)==4
                            Pcut4 = [Pcut4;xi1(j,:)];
                        elseif statepoint(j)==-1
                            Pin = [Pin;xi1(j,:)];
                        elseif statepoint(j)==1
                            Pout = [Pout;xi1(j,:)];
                        end
                    end
                end
                if ~isempty(Pcut2)
                    plot(POINT(Pcut2),'d','MarkerSize',10,'Color','b','LineWidth',2)
                end
                if ~isempty(Pcut3)
                    plot(POINT(Pcut3),'x','MarkerSize',10,'Color','r','LineWidth',2)
                end
                if ~isempty(Pcut4)
                    plot(POINT(Pcut4),'h','MarkerSize',10,'Color','c','LineWidth',2)
                end
                if ~isempty(Pin)
                    plot(POINT(Pin),'square','MarkerSize',10,'Color','y','LineWidth',2)
                end
                if ~isempty(Pout)
                    plot(POINT(Pout),'o','MarkerSize',10,'Color','g','LineWidth',2);
                end
                hold off
            end
        case 3
            
            hold on
            nbelem = length(H.e);
            stateH = getstate(H);
            
            for i =1:nbelem
                elem = H.e{i};
                xi1=calc_sommets(calc_xi_ani(elem,3),3);
                xi = [xi1(1,:);xi1(2,:);xi1(4,:);xi1(3,:);xi1(5,:);xi1(6,:);xi1(8,:);xi1(7,:);];
                faces = [ 1 2 6 5;2 3 7 6; 3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
                
                if state==1
                    col = stateH(i);
                    options = {'FaceColor' 'flat' 'FaceVertexCData' col};
                    options = [options , varargin];
                else
                    options = {'FaceColor'  'none'};
                    options = [options , varargin];
                end
                
                patch('Faces',faces,'Vertices',xi,options{:})
                
            end
            view(3); axis square
            
    end
else
    
    matdecoup = getcharin('matdecoup',varargin);
    varargin = delonlycharin('matdecoup',varargin);
    stodim = getM(H);
    vectdecoup = setdiff(1:stodim,matdecoup(:,1));
    
    
    hold on
    view(3); axis square
    xlim([0 1])
    ylim([0 1])
    nbelem = length(H.e);
    stateH = getstate(H);
    
    for i =1:nbelem
        elem = H.e{i};
        maxorder = elem.order;
        test = zeros(size(matdecoup(:,2)));
        for j=1:(stodim-3)
            decoup = 0:1/2^(maxorder(matdecoup(j,1))):1;
            test(j) = decoup(find(decoup>=matdecoup(j,2) & decoup<=matdecoup(j,2)+1/2^(maxorder(matdecoup(:,1)))));
        end
        
        I = calc_multi_indices(3,1,2);
        I=I(:,1:end-1);
        xi = calc_xi_ani(elem,stodim);
        %testxi1 = xi(1,matdecoup(:,1));
        testxi2 = xi(2,matdecoup(:,1));
        %keyboard
        if  all(test==testxi2)
            for j=1:3
                for k=1:8
                    if I(k,j)==0
                        I(k,j) = xi(1,vectdecoup(j));
                    elseif I(k,j)==1
                        I(k,j) = xi(2,vectdecoup(j));
                    end
                end
            end
            vertex = [I(1,:);I(2,:);I(4,:);I(3,:);I(5,:);I(6,:);I(8,:);I(7,:)];
            faces = [ 1 2 6 5;2 3 7 6; 3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
            if state==1
                col = stateH(i);
                options = {'FaceColor' 'flat' 'FaceVertexCData' col};
                
            else
                options = {'FaceColor'  'none'};
                
            end
            patch('Faces',faces,'Vertices',vertex,options{:})
        end
    end
    view(3); axis square
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

