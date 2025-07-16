function mov = evolLocalSolutionReference(patches,wt_ref,varargin)
% function mov = evolLocalSolutionReference(patches,wt_ref,varargin)
% Display evolution of reference local solution wt_ref
% patches: Patches or Patch
% wt_ref: TIMEMATRIX of reference local solution w
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
        figure('Name',['Reference local solution sig_w_ref_' num2str(i) ' over patches #' num2str([numbers{:}])])
        % set(gcf,'Name',['Reference local solution sig_w_ref_' num2str(i) ' over patches #' num2str([numbers{:}])])
    elseif ischarin('epsilon',varargin)
        i = getcharin('epsilon',varargin);
        figure('Name',['Reference local solution eps_w_ref_' num2str(i) ' over patches #' num2str([numbers{:}])])
        % set(gcf,'Name',['Reference local solution eps_w_ref_' num2str(i) ' over patches #' num2str([numbers{:}])])
    elseif ischarin('energyint',varargin)
        figure('Name',['Reference local solution H_w_ref over patches #' num2str([numbers{:}])])
        % set(gcf,'Name',['Reference local solution H_w_ref over patches #' num2str([numbers{:}])])
    elseif ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Reference local solution w_ref_' num2str(i) ' over patches #' num2str([numbers{:}])])
        % set(gcf,'Name',['Reference local solution w_ref_' num2str(i) ' over patches #' num2str([numbers{:}])])
    elseif ischarin('rotation',varargin)
        i = getcharin('rotation',varargin);
        figure('Name',['Reference local solution rot_w_ref_' num2str(i) ' over patches #' num2str([numbers{:}])])
        % set(gcf,'Name',['Reference local solution rot_w_ref_' num2str(i) ' over patches #' num2str([numbers{:}])])
    else
        figure('Name',['Reference local solution w_ref over patches #' num2str([numbers{:}])])
        % set(gcf,'Name',['Reference local solution w_ref over patches #' num2str([numbers{:}])])
    end
    clf
    set(gcf,'Color','w')
    
    n = numel(patches);
    S = cellfun(@(patch) patch.S,patches.patches,'UniformOutput',false);
    for k=1:n
        patch = patches.patches{k};
        wt_ref{patch.number} = setevolparam(wt_ref{patch.number},'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,...
            'view',p.Results.view,'camup',p.Results.camup,'campos',p.Results.campos,'FontSize',p.Results.FontSize);
    end
    T = gettimemodel(wt_ref{1});
    frame = evol_sol(T,wt_ref,S,'rescale',p.Results.rescale,varargin{:}); % save the frames
elseif isa(patches,'Patch')
    patch = patches;
    if ischarin('sigma',varargin)
        i = getcharin('sigma',varargin);
        figure('Name',['Reference local solution sig_w_ref_' num2str(i) ' over patch #' num2str(patch.number)])
        % set(gcf,'Name',['Reference local solution sig_w_ref_' num2str(i) ' over patch #' num2str(patch.number)])
    elseif ischarin('epsilon',varargin)
        i = getcharin('epsilon',varargin);
        figure('Name',['Reference local solution eps_w_ref_' num2str(i) ' over patch #' num2str(patch.number)])
        % set(gcf,'Name',['Reference local solution eps_w_ref_' num2str(i) ' over patch #' num2str(patch.number)])
    elseif ischarin('energyint',varargin)
        figure('Name',['Reference local solution H_w_ref over patch #' num2str(patch.number)])
        % set(gcf,'Name',['Reference local solution H_w_ref over patch #' num2str(patch.number)])
    elseif ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Reference local solution w_ref_' num2str(i) ' over patch #' num2str(patch.number)])
        % set(gcf,'Name',['Reference local solution w_ref_' num2str(i) ' over patch #' num2str(patch.number)])
    elseif ischarin('rotation',varargin)
        i = getcharin('rotation',varargin);
        figure('Name',['Reference local solution rot_w_ref_' num2str(i) ' over patch #' num2str(patch.number)])
        % set(gcf,'Name',['Reference local solution rot_w_ref_' num2str(i) ' over patch #' num2str(patch.number)])
    else
        figure('Name',['Reference local solution w_ref over patch #' num2str(patch.number)])
        % set(gcf,'Name',['Reference local solution w_ref over patch #' num2str(patch.number)])
    end
    clf
    set(gcf,'Color','w')
    
    wt_ref{patch.number} = setevolparam(wt_ref{patch.number},'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,...
        'view',p.Results.view,'camup',p.Results.camup,'campos',p.Results.campos,'FontSize',p.Results.FontSize);
    frame = evol_sol(wt_ref{patch.number},patch.S,'rescale',p.Results.rescale,varargin{:}); % save the frames
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
