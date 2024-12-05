function mov = evolClosedSobolIndicesLocalSolution(patches,T,wt,alpha,varargin)
% function mov = evolClosedSobolIndicesLocalSolution(patches,T,wt,alpha,varargin)
% Display evolution of the Closed Sobol indices of local solution wt
% associated with the group of variables alpha in {1,..,d} with d = ndims(wt)
% patches: Patches or Patch
% T: TIMEMODEL
% wt: FunctionalBasis of local solution w
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
addParameter(p,'FontSize',16,@isscalar);
addParameter(p,'filename','solution',@ischar);
addParameter(p,'pathname','./',@ischar);
addParameter(p,'formats',{'avi','mp4'},@(x) ischar(x) || iscell(x));
addParameter(p,'FrameRate',30,@isnumeric);
addParameter(p,'Quality',100,@isnumeric);
parse(p,varargin{:})

varargin = delcharin({'rescale','colorbar','colormap','view','camup','FontSize',...
    'filename','pathname','formats','FrameRate','Quality'},varargin);
if isa(p.Results.formats,'char')
    p.Results.formats = {p.Results.formats};
end

if isa(patches,'Patches')
    numbers = getnumber(patches);
    if ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Closed Sobol index of w_' num2str(i) ' over patches #' num2str([numbers{:}]) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Closed Sobol index of w_' num2str(i) ' over patches #' num2str([numbers{:}]) ' for random variables #' num2str(alpha)])
    else
        figure('Name',['Closed Sobol index of w over patches #' num2str([numbers{:}]) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Closed Sobol index of w over patches #' num2str([numbers{:}]) ' for random variables #' num2str(alpha)])
    end
    clf
    set(gcf,'color','w')
    
    patch = patches.patches;
    S = cellfun(@(patch) patch.S,patch,'UniformOutput',false);
    swt = cellfun(@(patch) TIMEMATRIX(reshape(SensitivityAnalysis.closedSobolIndices(wt{patch.number},alpha,ndims(wt{patch.number})),wt{patch.number}.sz),T),patch,'UniformOutput',false);
    T = setevolparam(T,'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,...
        'view',p.Results.view,'camup',p.Results.camup,'FontSize',p.Results.FontSize);
    frame = evol_sol(T,swt,S,'rescale',p.Results.rescale,varargin{:}); % save the frames
elseif isa(patches,'Patch')
    patch = patches;
    if ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Closed Sobol index of w_' num2str(i) ' over patch #' num2str(patch.number) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Closed Sobol index of w_' num2str(i) ' over patch #' num2str(patch.number) ' for random variables #' num2str(alpha)])
    else
        figure('Name',['Closed Sobol index of w over patch #' num2str(patch.number) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Closed Sobol index of w over patch #' num2str(patch.number) ' for random variables #' num2str(alpha)])
    end
    clf
    set(gcf,'color','w')
    
    sz = wt{patch.number}.sz;
    swt = SensitivityAnalysis.closedSobolIndices(wt{patch.number},alpha);
    swt = reshape(swt,sz);
    swt = TIMEMATRIX(swt,T);
    
    swt = setevolparam(swt,'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,...
        'view',p.Results.view,'camup',p.Results.camup,'FontSize',p.Results.FontSize);
    frame = evol_sol(swt,patch.S,'rescale',p.Results.rescale,varargin{:}); % save the frames
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
