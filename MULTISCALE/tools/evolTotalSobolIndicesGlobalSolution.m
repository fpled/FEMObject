function mov = evolTotalSobolIndicesGlobalSolution(glob,T,Ut,alpha,varargin)
% function mov = evolTotalSobolIndicesGlobalSolution(glob,T,Ut,alpha,varargin)
% Display evolution of the Total Sobol indices of global solution Ut
% associated with the group of variables alpha in {1,..,d} with d = ndims(Ut)
% glob: Global
% T: TIMEMODEL
% Ut: FunctionalBasisArray of global solution U
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

% Create a new figure
if ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Total Sobol index of U_' num2str(i) ' over fictitious domain for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Total Sobol index of U_' num2str(i) ' over fictitious domain for random variables #' num2str(alpha)])
else
    figure('Name',['Total Sobol index of U over fictitious domain for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Total Sobol index of U over fictitious domain for random variables #' num2str(alpha)])
end
clf
set(gcf,'color','w')

sz = Ut.sz;
d = ndims(Ut);
sUt = SensitivityAnalysis.totalSobolIndices(Ut,alpha,d);
sUt = reshape(sUt,sz);
sUt = TIMEMATRIX(sUt,T);
sUt = unfreevector(glob.S,sUt)-calc_init_dirichlet(glob.S)*one(T);

sUt = setevolparam(sUt,'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,...
    'view',p.Results.view,'camup',p.Results.camup,'campos',p.Results.campos,'FontSize',p.Results.FontSize);
frame = evol_sol(sUt,glob.S,'rescale',p.Results.rescale,varargin{:}); % save the frames

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
