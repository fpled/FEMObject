function [feModel,mesherTime] = buildFEModel(model,method)
% [feModel,mesherTime] = buildFEModel(model,method)

mesherClock = tic ;

if nargin == 1
    method = 1 ;
end

cellNum = getCellNum(model) ;
cellNb = getCellNb(model) ;
cellSize = getCellSize(model) ;

switch method
    case 1
        cellModel = getCellModel(model) ;
        feModel = MODEL() ;
        for j=1:cellNum(2)
            for i=1:cellNum(1)
                newModel = cellModel + ...
                    VECTOR([cellSize(1)*(i-1) ; cellSize(2)*(j-1)]) ;
                feModel = union(feModel,newModel) ;
                ifprint(model,sprintf('Period %i out of %i completed\n',...
                    i+(j-1)*cellNum(1),cellNb))
            end
        end
        
    case 2
        % Method works only for constant element size (including across patches)
        domainSize = cellSize.*cellNum ;
        D_S = DOMAIN(2,zeros(1,2),domainSize) ;
        elementSize = getElementSize(model) ;
        nbElem = round(domainSize./elementSize) ;
        nbElem = max(1,nbElem) ; % for safety
        feModel = mesh(D_S,nbElem(1),nbElem(2)) ;
        for j=1:cellNum(2)
            for i=1:cellNum(1)
                if ~(i==1 && j==1)
                    Dnew = DOMAIN(2,cellSize.*[(i-1) (j-1)],cellSize.*[i j]) ;
                    [~,~,numelem] = intersect(feModel,Dnew) ;
                    feModel = separateelemwithnum(feModel,numelem) ;
                    ifprint(model,sprintf('Period %i out of %i completed\n',...
                        i+(j-1)*cellNum(1),cellNb))
                end
            end
        end
end
feModel = createddlnode(feModel,DDL('u')) ;

mesherTime = toc(mesherClock) ;

end