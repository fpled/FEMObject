function [feModel,mesherTime] = buildFEModel(pb,method)
% [feModel,mesherTime] = buildFEModel(pb,method)

mesherClock = tic ;

if pb.coordTol~=getfemobjectoptions('tolerancepoint')
    warning('Problem''s tolerance on coordinate differs from femObject''s.')
end

if nargin == 1
    method = 1 ;
end

groupelemByCell = false ; % hard-coded switch

model = pb.model ;
cellNum = getCellNum(model) ;
cellNb = getCellNb(model) ;
cellSize = getCellSize(model) ;
patchNb = numel(pb.patch) ;
% List lower left cell of each patch
patchLLCellList = zeros(patchNb,getOrder(pb.model)-1) ;
for p = 1:patchNb
    patchLLCellList(p,:) = formatIndex(getOrder(pb.model),cellNum,...
        min(formatIndex(2,cellNum,pb.patchCells{p}))) ;
end

if method == 1
    cellModel = getCellModel(model) ;
    patchCellList = formatIndex(2,cellNum,allPatchCells(pb)) ;
    feModel = MODEL() ;
    for j=1:cellNum(2)
        for i=1:cellNum(1)
            ind = sub2ind(cellNum,i,j) ;
            if ~ismember(ind,patchCellList,'rows')
                newModel = cellModel ;
            elseif ismember(ind,patchLLCellList,'rows') 
                % Find patch model
                [~,repPatch] = ismember(ind,patchLLCellList,'rows') ;
                ifprint(pb.model.verbose,...
                    'Cell %i is lower left in patch %i. Building FE model...\n',...
                    ind,repPatch)
                % Build patch model
                newModel = buildFEModel(pb.patch(repPatch)) ;
                ifprint(pb.model.verbose,'FE model for patch %i built.\n',repPatch)
            else % already in a patch model
                ifprint(pb.model.verbose,'Skipped cell %i out of %i: part of a patch.\n',ind,cellNb)
                continue
            end
            newModel = newModel + ...
                VECTOR([cellSize(1)*(i-1) ; cellSize(2)*(j-1)]) ;
            feModel = union(feModel,newModel) ;
            ifprint(pb.model.verbose,sprintf('Cell %i out of %i completed\n',...
                ind,cellNb))
        end
    end
    
elseif method == 2 
    % Method works only for constant element size (including across patches)
    ifprint(pb.model.verbose,'Building global mesh of')
    domainSize = cellSize.*cellNum ;
    D_S = DOMAIN(2,zeros(1,2),domainSize) ;
    elementSize = getElementSize(model) ;
    nbElem = round(domainSize./elementSize) ;
    nbElem = max(1,nbElem) ; % for safety
    ifprint(pb.model.verbose,' %i by %i elements...',nbElem(1),nbElem(2))
    feModel = mesh(D_S,nbElem(1),nbElem(2)) ;
    ifprint(pb.model.verbose,' Done.\n')
    if groupelemByCell
        ifprint(pb.model.verbose,'Will group elements by cell.\n')
        for j=1:cellNum(2)
            for i=1:cellNum(1)
                if ~(i==1 && j==1)
                    Dnew = DOMAIN(2,cellSize.*[(i-1) (j-1)],cellSize.*[i j]) ;
                    [~,~,numelem] = intersect(feModel,Dnew) ;
                    feModel = separateelemwithnum(feModel,numelem) ;
                    ifprint(model.verbose,sprintf('Cell %i out of %i completed\n',...
                        i+(j-1)*cellNum(1),cellNb))
                end
            end
        end
    end
    for p = 1:patchNb
        ifprint(pb.model.verbose,...
            'Local mesh substitution for patch %i of %i.\n',p,patchNb)
        % Build patch model
        ifprint(pb.model.verbose,'Building patch model...\n')
        patchModel = buildFEModel(pb.patch(p),1) ;
        % Get coordinates of lower left corner
        ll = patchLLCellList(p,:) ;
        ll = (formatIndex(3,cellNum,ll)-1)*diag(getCellSize(pb.model)) ;
        patchModel = patchModel + VECTOR(ll(:)) ;
        ifprint(pb.model.verbose,' Done.\n')
        % Substitute
        ifprint(pb.model.verbose,'Substituting mesh on patch...')
        patchDomain = DOMAIN(2,ll,ll+getDomainSize(pb.patch(p).model)) ;
        feModel = setdiff(feModel,patchDomain,'strict',0) ;
        feModel = union(feModel,patchModel) ;
        ifprint(pb.model.verbose,' Done.\n')
    end    
else
    error('Method unknown')
end
feModel = createddlnode(feModel,DDL('u')) ;

mesherTime = toc(mesherClock) ;

end