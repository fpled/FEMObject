function mov = evolStdLocalSolutionReference(patches,T,wt_ref,varargin)
% function mov = evolStdLocalSolutionReference(patches,T,wt_ref,varargin)
% Display evolution of the standard deviation of reference local solution wt_ref
% patches: Patches or Patch
% T: TIMEMODEL
% wt_ref: FunctionalBasis of reference local solution w
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
    if ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Standard deviation of w_ref_' num2str(i) ' over patches #' num2str([numbers{:}])])
        % set(gcf,'Name',['Standard deviation of w_ref_' num2str(i) ' over patches #' num2str([numbers{:}])])
    else
        figure('Name',['Standard deviation of w_ref over patches #' num2str([numbers{:}])])
        % set(gcf,'Name',['Standard deviation of w_ref over patches #' num2str([numbers{:}])])
    end
    clf
    set(gcf,'color','w')
    
    patch = patches.patches;
    S = cellfun(@(patch) patch.S,patch,'UniformOutput',false);
    swt_ref = cellfun(@(patch) TIMEMATRIX(reshape(variance(wt_ref{patch.number}),wt_ref{patch.number}.sz),T),patch,'UniformOutput',false);
    T = setevolparam(T,'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,...
        'view',p.Results.view,'camup',p.Results.camup,'campos',p.Results.campos,'FontSize',p.Results.FontSize);
    frame = evol_sol(T,swt_ref,S,'rescale',p.Results.rescale,varargin{:}); % save the frames
elseif isa(patches,'Patch')
    patch = patches;
    if ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Standard deviation of w_ref_' num2str(i) ' over patch #' num2str(patch.number)])
        % set(gcf,'Name',['Standard deviation of w_ref_' num2str(i) ' over patch #' num2str(patch.number)])
    else
        figure('Name',['Standard deviation of w_ref over patch #' num2str(patch.number)])
        % set(gcf,'Name',['Standard deviation of w_ref over patch #' num2str(patch.number)])
    end
    clf
    set(gcf,'color','w')
    
    sz = wt_ref{patch.number}.sz;
    swt_ref = variance(wt_ref{patch.number});
    swt_ref = reshape(swt_ref,sz);
    swt_ref = TIMEMATRIX(swt_ref,T);
    
    swt_ref = setevolparam(swt_ref,'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,...
        'view',p.Results.view,'camup',p.Results.camup,'campos',p.Results.campos,'FontSize',p.Results.FontSize);
    frame = evol_sol(swt_ref,patch.S,'rescale',p.Results.rescale,varargin{:}); % save the frames
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
