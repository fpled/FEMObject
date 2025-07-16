function mov = evolSobolIndicesGlobalLocalSolution(glob,patches,interfaces,T,Ut,wt,alpha,varargin)
% function mov = evolSobolIndicesGlobalLocalSolution(glob,patches,interfaces,T,Ut,wt,alpha,varargin)
% Display evolution of the Sobol indices of global solution Ut and local solution wt
% associated with the group of variables alpha in {1,..,d} with d = min(ndims(Ut),ndims(wt))
% glob: Global
% patches: Patches
% interfaces: Interfaces
% T: TIMEMODEL
% Ut: TIMEMATRIX of global solution U
% wt: TIMEMATRIX of local solution w
% alpha: 1-by-s array of integers or 1-by-d logical
% - if alpha is an array of integers, indices with respect
% to variables alpha
% - if alpha is logical, indices with respect
% to variables find(alpha)
% mov: movie

p = ImprovedInputParser;
addParameter(p,'rescale',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'view',[],@isnumeric);
addParameter(p,'camup','auto',@(x) isnumeric(x) || ischar(x));
addParameter(p,'campos','auto',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
addParameter(p,'filename','solution',@ischar);
addParameter(p,'pathname','./',@ischar);
addParameter(p,'formats',{'avi','mp4'},@(x) ischar(x) || iscell(x));
addParameter(p,'FrameRate',30,@isnumeric);
addParameter(p,'Quality',100,@isnumeric);
parse(p,varargin{:})

varargin = delcharin({'rescale','colorbar','colormap','view','camup','campos','FontSize',...
    'filename','pathname','formats','FrameRate','Quality'},varargin);
if isa(p.Results.formats,'char')
    p.Results.formats = {p.Results.formats};
end
n = numel(patches);

if ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Sobol index of U_' num2str(i) ' and w_' num2str(i) ' over domain and patches for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Sobol index of U_' num2str(i) ' and w_' num2str(i) ' over domain and patches for random variables #' num2str(alpha)])
else
    figure('Name',['Sobol index of U and w over domain and patches for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Sobol index of U and w over domain and patches for random variables #' num2str(alpha)])
end
clf
set(gcf,'Color','w')

sz = Ut.sz;
dUt = ndims(Ut);
sUt = SensitivityAnalysis.sobolIndices(Ut,alpha,dUt);
sUt = reshape(sUt,sz);
sUt = TIMEMATRIX(sUt,T);
sUt = unfreevector(glob.S,sUt)-calc_init_dirichlet(glob.S)*one(T);

swt = cellfun(@(patch) TIMEMATRIX(reshape(SensitivityAnalysis.sobolIndices(wt{patch.number},alpha,ndims(wt{patch.number})),wt{patch.number}.sz),T),patches.patches,'UniformOutput',false);

T = setevolparam(T,'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,...
    'view',p.Results.view,'camup',p.Results.camup,'campos',p.Results.campos,'FontSize',p.Results.FontSize);
if p.Results.view==3
    T = setevolparam(T,'axisimage',true,'axissquare',false);
end
[t,rep]=gettevol(T);

if p.Results.rescale
    
    qtmax = cell(1,1+n);
    qtmin = cell(1,1+n);
    if ischarin('sigma',varargin)
        q = calc_sigma(glob.S,sUt,varargin{2:end});
        k = getcharin('sigma',varargin);
        q = sigmacompo(q,k,glob.S);
    elseif ischarin('epsilon',varargin)
        q = calc_epsilon(glob.S,sUt,varargin{2:end});
        k = getcharin('epsilon',varargin);
        q = sigmacompo(q,k,glob.S);
    elseif ischarin('energyint',varargin)
        q = calc_energyint(glob.S,sUt,varargin{2:end});
    elseif ischarin('displ',varargin)
        q = unfreevector(glob.S,sUt);
        k = getcharin('displ',varargin);
        switch k
            case 1
                numddl = findddl(glob.S,'UX');
            case 2
                numddl = findddl(glob.S,'UY');
            case 3
                numddl = findddl(glob.S,'UZ');
        end
        q = getcompo(q,numddl);
    elseif ischarin('rotation',varargin)
        q = unfreevector(glob.S,sUt);
        k = getcharin('rotation',varargin);
        switch k
            case 1
                numddl = findddl(glob.S,'RX');
            case 2
                numddl = findddl(glob.S,'RY');
            case 3
                numddl = findddl(glob.S,'RZ');
        end
        q = getcompo(q,numddl);
    else
        q = sUt;
    end
    
    if isa(q,'TIMEMATRIX') && ~israndom(q)
        qtmax{1} = max(max(q));
        qtmin{1} = min(min(q));
    elseif isa(q,'TIMERADIALMATRIX')
        qtmax{1}=-Inf;
        qtmin{1}=Inf;
        for kk=1:max(1,floor(length(rep)/5)):length(rep)
            qtmax{1} = max(qtmax{1},max(max(getmatrixatstep(q,rep(kk)))));
            qtmin{1} = min(qtmin{1},min(min(getmatrixatstep(q,rep(kk)))));
        end
    end
    
    for j=1:n
        patch = patches.patches{j};
        if ischarin('sigma',varargin)
            q = calc_sigma(patch.S,swt{patch.number},varargin{2:end});
            k = getcharin('sigma',varargin);
            q = sigmacompo(q,k,patch.S);
        elseif ischarin('epsilon',varargin)
            q = calc_epsilon(patch.S,swt{patch.number},varargin{2:end});
            k = getcharin('epsilon',varargin);
            q = sigmacompo(q,k,patch.S);
        elseif ischarin('energyint',varargin)
            q = calc_energyint(patch.S,swt{patch.number},varargin{2:end});
        elseif ischarin('displ',varargin)
            q = unfreevector(patch.S,swt{patch.number});
            k = getcharin('displ',varargin);
            switch k
                case 1
                    numddl = findddl(patch.S,'UX');
                case 2
                    numddl = findddl(patch.S,'UY');
                case 3
                    numddl = findddl(patch.S,'UZ');
            end
            q = getcompo(q,numddl);
        elseif ischarin('rotation',varargin)
            q = unfreevector(patch.S,swt{patch.number});
            k = getcharin('rotation',varargin);
            switch k
                case 1
                    numddl = findddl(patch.S,'RX');
                case 2
                    numddl = findddl(patch.S,'RY');
                case 3
                    numddl = findddl(patch.S,'RZ');
            end
            q = getcompo(q,numddl);
        else
            q = swt{patch.number};
        end
        
        if isa(q,'TIMEMATRIX') && ~israndom(q)
            qtmax{j+1} = max(max(q));
            qtmin{j+1} = min(min(q));
        elseif isa(q,'TIMERADIALMATRIX')
            qtmax{j+1}=-Inf;
            qtmin{j+1}=Inf;
            for kk=1:max(1,floor(length(rep)/5)):length(rep)
                qtmax{j+1} = max(qtmax{j+1},max(max(getmatrixatstep(q,rep(kk)))));
                qtmin{j+1} = min(qtmin{j+1},min(min(getmatrixatstep(q,rep(kk)))));
            end
        end
    end
    qtmax = max([qtmax{:}]);
    qtmin = min([qtmin{:}]);
    
    T = setevolparam(T,'caxis',[qtmin,qtmax]);
    if strfind(p.Results.rescale,'z')
        T = setevolparam(T,'zlim',[qtmin,qtmax]);
    elseif strfind(p.Results.rescale,'y')
        T = setevolparam(T,'ylim',[qtmin,qtmax]);
    end
    varargin = delcharin('rescale',varargin);
    
end

plotiter = getevolparam(T,'plotiter');
plottime = getevolparam(T,'plottime');
plotstep = getevolparam(T,'plotstep');
colorba = getevolparam(T,'colorbar');
compt = 0;
for i=1:plotstep:length(t)
    compt=compt+1;
    execute(T.evolparam,'before');
    
    h1 = subplot(2,1,1);
    for k=1:n
        patch = patches.patches{k};
        interface = interfaces.interfaces{k};
        plot_sol(patch.S,getmatrixatstep(swt{patch.number},rep(i)),varargin{:});
        if ~ischarin('sigma',varargin) && ~ischarin('epsilon',varargin) && ~ischarin('energyint',varargin)
            plot_sol(interface.S,interface.P_patch*getmatrixatstep(swt{patch.number},rep(i)),'FaceColor','none','EdgeColor','k',varargin{:});
        end
    end
    if colorba
        c1 = colorbar;
        posc1 = get(c1,'Position');
        posh1 = get(h1,'Position');
        colorbar('off');
    end
    ax1 = axis;
    cax1 = caxis;
    if colorba
        set(h1,'Position',posh1);
    end
    
    TT = setevolparam(T,'pause',false);
    TT = setevolparam(TT,'colorbar',false);
    execute(TT.evolparam,'after')
    
    h2 = subplot(2,1,2);
    plot_sol(glob.S,getmatrixatstep(sUt,rep(i)),varargin{:});
    for k=1:n
        patch = patches.patches{k};
        interface = interfaces.interfaces{k};
        if ~ischarin('sigma',varargin) && ~ischarin('epsilon',varargin) && ~ischarin('energyint',varargin)
            plot_sol(interface.S,interface.P_patch*getmatrixatstep(swt{patch.number},rep(i)),'FaceColor','none','EdgeColor','k',varargin{:});
        end
    end
    if colorba
        c2 = colorbar;
        posc2 = get(c2,'Position');
        posh2 = get(h2,'Position');
        posc2(4) = posc1(2)+posc1(4)-posc2(2);
        set(c2,'Position',posc2);
    end
    ax2 = axis;
    cax2 = caxis;
    if colorba
        set(h2,'Position',posh2);
    end
    
    ax = ax1;
    ax(1:2:end) = min(ax1(1:2:end),ax2(1:2:end));
    ax(2:2:end) = max(ax1(2:2:end),ax2(2:2:end));
    % if getevolparam(T,'view')==3
    %     ax(5:6) = [-eps eps];
    % end
    axis([h1 h2],ax)
    
    if ~p.Results.rescale
        cax(1) = min(cax1(1),cax2(1));
        cax(2) = max(cax1(2),cax2(2));
        caxis(h1,cax)
        caxis(h2,cax)
    end
    
    pos1 = get(h1,'Position');
    pos2 = get(h2,'Position');
    trans = pos1(2)-(pos2(2)+pos2(4));
    pos1(2) = pos1(2)-trans;
    set(h1,'Position',pos1);
    
    if colorba
        posc2(4) = posc2(4)-trans;
        set(c2,'Position',posc2);
    end
    
    if plotiter
        title(h1,['iter ' num2str(i,'%d')],'FontSize',p.Results.FontSize)
    elseif plottime
        title(h1,['time ' num2str(t(i),'%.2f') ' s'],'FontSize',p.Results.FontSize)
    end
    
    T = setevolparam(T,'colorbar',false);
    execute(T.evolparam,'after')
    
    frame(compt) = getframe(gcf);
    
end

% Create movie file
mov = cell(1,length(p.Results.formats));
for i=1:length(p.Results.formats)
    if strcmp(p.Results.formats{i},'avi')
        mov{i} = VideoWriter(fullfile(p.Results.pathname,p.Results.filename));
    elseif strcmp(p.Results.formats{i},'mp4')
        mov{i} = VideoWriter(fullfile(p.Results.pathname,p.Results.filename),'MPEG-4');
    elseif strcmp(p.Results.formats{i},'mj2')
        mov{i} = VideoWriter(fullfile(p.Results.pathname,p.Results.filename),'Motion JPEG 2000');
    end
    mov{i}.FrameRate = p.Results.FrameRate;
    mov{i}.Quality = p.Results.Quality;
    open(mov{i});
    writeVideo(mov{i},frame); % add the frames to the movie
    close(mov{i});
end


function se = sigmacompo(se,ksigma,S)

if isa(ksigma,'char') && strcmp(ksigma,'mises')
    switch getindim(S)
        case 1
            se = getcompo(se,1);
        case 2
            se1 = getcompo(se,1);
            se2 = getcompo(se,2);
            se3 = getcompo(se,3);
            tracese = 1/3*(se1+se2);
            se1 = se1 - tracese;
            se2 = se2 - tracese;
            T = gettimemodel(se);
            se = sqrt(3/2*(getvalue(se1).*getvalue(se1) + getvalue(se2).*getvalue(se2) + 2*getvalue(se3).*getvalue(se3)));
            se = TIMEMATRIX(se,T);
        case 3
            se1 = getcompo(se,1);
            se2 = getcompo(se,2);
            se3 = getcompo(se,3);
            se4 = getcompo(se,4);
            se5 = getcompo(se,5);
            se6 = getcompo(se,6);
            tracese = 1/3*(se1+se2+se3);
            se1 = se1 - tracese;
            se2 = se2 - tracese;
            se3 = se3 - tracese;
            T = gettimemodel(se);
            se = sqrt(3/2*(getvalue(se1).*getvalue(se1) + getvalue(se2).*getvalue(se2) + getvalue(se3).*getvalue(se3)...
                + 2*(getvalue(se4).*getvalue(se4) + getvalue(se5).*getvalue(se5) + getvalue(se6).*getvalue(se6))));
            se = TIMEMATRIX(se,T);
    end
else
    se = getcompo(se,ksigma);
end

return
