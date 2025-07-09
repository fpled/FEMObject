function mov = evolLocalSolution(patches,wt,varargin)
% function mov = evolLocalSolution(patches,wt,varargin)
% Display evolution of local solution wt
% patches: Patches or Patch
% wt: TIMEMATRIX of local solution w
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

if isa(patches,'Patches')
    numbers = getnumber(patches);
    if ischarin('sigma',varargin)
        i = getcharin('sigma',varargin);
        figure('Name',['Local solution sig_w_' num2str(i) ' over patches #' num2str([numbers{:}])])
        % set(gcf,'Name',['Local solution sig_w_' num2str(i) ' over patches #' num2str([numbers{:}])])
    elseif ischarin('epsilon',varargin)
        i = getcharin('epsilon',varargin);
        figure('Name',['Local solution eps_w_' num2str(i) ' over patches #' num2str([numbers{:}])])
        % set(gcf,'Name',['Local solution eps_w_' num2str(i) ' over patches #' num2str([numbers{:}])])
    elseif ischarin('energyint',varargin)
        figure('Name',['Local solution H_w over patches #' num2str([numbers{:}])])
        % set(gcf,'Name',['Local solution H_w over patches #' num2str([numbers{:}])])
    elseif ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Local solution w_' num2str(i) ' over patches #' num2str([numbers{:}])])
        % set(gcf,'Name',['Local solution w_' num2str(i) ' over patches #' num2str([numbers{:}])])
    elseif ischarin('rotation',varargin)
        i = getcharin('rotation',varargin);
        figure('Name',['Local solution rot_w_' num2str(i) ' over patches #' num2str([numbers{:}])])
        % set(gcf,'Name',['Local solution rot_w_' num2str(i) ' over patches #' num2str([numbers{:}])])
    else
        figure('Name',['Local solution w over patches #' num2str([numbers{:}])])
        % set(gcf,'Name',['Local solution w over patches #' num2str([numbers{:}])])
    end
    clf
    set(gcf,'color','w')
    
    n = numel(patches);
    S = cellfun(@(patch) patch.S,patches.patches,'UniformOutput',false);
    for k=1:n
        patch = patches.patches{k};
        wt{patch.number} = setevolparam(wt{patch.number},'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,...
            'view',p.Results.view,'camup',p.Results.camup,'campos',p.Results.campos,'FontSize',p.Results.FontSize);
    end
    T = gettimemodel(wt{1});
    frame = evol_sol(T,wt,S,'rescale',p.Results.rescale,varargin{:}); % save the frames
elseif isa(patches,'Patch')
    patch = patches;
    if ischarin('sigma',varargin)
        i = getcharin('sigma',varargin);
        figure('Name',['Local solution sig_w_' num2str(i) ' over patch #' num2str(patch.number)])
        % set(gcf,'Name',['Local solution sig_w_' num2str(i) ' over patch #' num2str(patch.number)])
    elseif ischarin('epsilon',varargin)
        i = getcharin('epsilon',varargin);
        figure('Name',['Local solution eps_w_' num2str(i) ' over patch #' num2str(patch.number)])
        % set(gcf,'Name',['Local solution eps_w_' num2str(i) ' over patch #' num2str(patch.number)])
    elseif ischarin('energyint',varargin)
        figure('Name',['Local solution H_w over patch #' num2str(patch.number)])
        % set(gcf,'Name',['Local solution H_w over patch #' num2str(patch.number)])
    elseif ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Local solution w_' num2str(i) ' over patch #' num2str(patch.number)])
        % set(gcf,'Name',['Local solution w_' num2str(i) ' over patch #' num2str(patch.number)])
    elseif ischarin('rotation',varargin)
        i = getcharin('rotation',varargin);
        figure('Name',['Local solution rot_w_' num2str(i) ' over patch #' num2str(patch.number)])
        % set(gcf,'Name',['Local solution rot_w_' num2str(i) ' over patch #' num2str(patch.number)])
    else
        figure('Name',['Local solution w over patch #' num2str(patch.number)])
        % set(gcf,'Name',['Local solution w over patch #' num2str(patch.number)])
    end
    clf
    set(gcf,'color','w')
    
    wt{patch.number} = setevolparam(wt{patch.number},'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,'FontSize',p.Results.FontSize);
    frame = evol_sol(wt{patch.number},patch.S,'rescale',p.Results.rescale,varargin{:}); % save the frames
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

end
