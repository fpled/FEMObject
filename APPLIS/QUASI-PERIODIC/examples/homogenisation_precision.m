% Parameters

caseChoice = 1 ; % 1 unidirectional | 2 bidirectional
maxIter = 20 ;
stagnationTol = 0 ;
approxTol = [1e-1 1e-3 1e-5] ;
mesoSize = 400 ;
contrastK = 100 ;
proba = .1 ;

%% Pre-processing

switch caseChoice
    case 1
        model = QPModel('order',2,...
            'cellNum',[mesoSize(1) 1],...
            'cellSize',[1 mesoSize(1)],...
            'elementSize',[1 mesoSize(1)]/20,...
            'tolSVD',1e-9,...
            'verbose',false) ;
        patterns = struct('name',{'uniform','bar'},...
            'value',{1 contrastK(1)-1},...
            'size',{[] [.5 1]} ,...
            'center',{[] [.5 .5]},...
            'offset',{[] []}) ;
        patternsTable = [1 1 ; 1 0] ;
                
    case 2
        model = QPModel('order',2,...
            'cellNum',floor(sqrt(mesoSize(1)*[1 1])),...
            'cellSize',[1 1],...
            'elementSize',[.05 .05],...
            'tolSVD',1e-9,...
            'verbose',false) ;
        patterns = struct('name',{'uniform','rectangle'},...
            'value',{1 contrastK(1)-1},...
            'size',{[] sqrt([.5 .5])} ,...
            'center',{[] [.5 .5]},...
            'offset',{[] []}) ;
        patternsTable = [1 1 ; 1 0] ;
end

KAss = QPConductivityAssembler('model',model,...
    'patterns',patterns,...
    'patternsTable',patternsTable,...
    'probability',[1-proba(1) proba(1)]) ;

diffAss = QPDiffusionAssembler('conductivityAssembler',KAss,...
    'BC','PBC',...
    'source','corrector1',...
    'useStabilisationWeights',true,...% required if caseChoice==1
    'useAverageWeights',true) ;       % for low tolerances (why ?)

pb = QPDiffusionProblem('operatorAssembler',diffAss,'tolerance',approxTol(3)) ;

% Build randomly maxIter conductivity fields
conductivities = cell(maxIter,1) ;
cellCoord = getCellCoord(model) ;
microPatterns = drawCellPattern(cellCoord,patterns)*patternsTable ;
for i = 1:maxIter
    dist = dealMultinomial([1-proba(1) proba(1)],getCellNb(model)) ;
    conductivities{i} = distributeMicro(model,dist,microPatterns) ;
end

% QP approximation
nParam = numel(approxTol) ;
outQP = cell(1,nParam) ;
for j = 1:nParam
    pb = setTolerance(pb,approxTol(j)) ;
    [KhomoQP,outQP{j}] = homogenize(pb,1,maxIter,stagnationTol,conductivities) ;
end
iter = max(cellfun(@(s)numel(s.times),outQP)) ;
outQP = cat(1,outQP{:}) ;

% FE computations
[KhomoFE,outFE] = homogenize(pb,2,iter,0,conductivities(1:iter)) ;

%% Post-processing

% matlab2tikz additional arguments
m2tikzArgin = {'showInfo',false,'showwarnings',false,'noSize',true,...
    'extraAxisOptions',{'ymajorgrids=true','grid style ={dotted}'}} ;

% Line specifications and legend
specP = {':','--','-'} ;
leg = cell(1+nParam,1) ;
spec = leg ;
leg{end} = 'FEM' ;
spec{end} =  [specP{end},'xr'] ;
for j = 1:nParam
    spec{j} = [specP{j},'ob'] ;
    leg{j} = sprintf('LR (precision %g)',approxTol(j)) ;
end

% Stagnation and time
figure
subplot(2,1,1)
hold on
for j = 1:nParam
    plot(2:numel(outQP(j).times),outQP(j).stagnations,spec{j})
end
plot(2:iter,outFE.stagnations,spec{end})
hold off
ylabel('Stagnation')
legend(leg(:))
subplot(2,1,2)
hold on
for j = 1:nParam
    plot(1:numel(outQP(j).times),log10(outQP(j).times),spec{j})
end
plot(1:iter,log10(outFE.times),spec{end})
hold off
ylabel('Time (log(s))')
legend(leg(:))
cleanfigure('pruneText',false); 
tikzName = sprintf('precision%iD_stag-time.tikz',caseChoice) ;
matlab2tikz(tikzName,m2tikzArgin{:})


% Average time with respect to precision
figure
hold on
avgTimes = zeros(nParam,1) ;
for j = 1:nParam
    avgTimes(j,1) = mean(outQP(j).times) ;
end
plot(-log10(approxTol)',avgTimes(:,1),'-xr')
plot(-log10(approxTol)',repmat(mean(outFE.times),nParam,1),'-ob')
hold off
ylabel('Average time (s)')
xlabel('-log(tolerance)')
legend({'FEM','LR'})
cleanfigure('pruneText',false); 
tikzName = sprintf('precision%iD_avgTime.tikz',caseChoice) ;
matlab2tikz(tikzName,m2tikzArgin{:})


% Homogenized coefficients convergence (V1)
figure
if caseChoice == 1
    subplot(2,1,1)
    hold on
    for j = 1:nParam
        plot(1:numel(outQP(j).times),cellfun(@(c)c(1,1),...
            outQP(j).effectiveConductivities),spec{j})
    end
    plot(1:iter,cellfun(@(c)c(1,1),...
        outFE.effectiveConductivities),spec{end})
    hold off
    ylabel('K_{11}')
    legend(leg(:))
    subplot(2,1,2)
end
hold on
for j = 1:nParam
    plot(1:numel(outQP(j).times),cellfun(@(c)c(2,2),...
        outQP(j).effectiveConductivities),spec{j})
end
plot(1:iter,cellfun(@(c)c(2,2),outFE.effectiveConductivities),spec{end})
hold off
ylabel('K_{22}')
legend(leg(:))
cleanfigure('pruneText',false); 
tikzName = sprintf('precision%iD_K.tikz',caseChoice) ;
matlab2tikz(tikzName,m2tikzArgin{:})