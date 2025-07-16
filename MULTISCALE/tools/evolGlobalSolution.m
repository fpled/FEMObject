function mov = evolGlobalSolution(glob,Ut,varargin)
% function mov = evolGlobalSolution(glob,Ut,varargin)
% Display evolution of global solution Ut
% glob: Global
% Ut: TIMEMATRIX of global solution U
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

if ischarin('sigma',varargin)
    i = getcharin('sigma',varargin);
    figure('Name',['Global solution sig_U_' num2str(i) ' over fictitious domain'])
    % set(gcf,'Name',['Global solution sig_U_' num2str(i) ' over fictitious domain'])
elseif ischarin('epsilon',varargin)
    i = getcharin('epsilon',varargin);
    figure('Name',['Global solution eps_U_' num2str(i) ' over fictitious domain'])
    % set(gcf,'Name',['Global solution eps_U_' num2str(i) ' over fictitious domain'])
elseif ischarin('energyint',varargin)
    figure('Name',['Global solution H_U over fictitious domain'])
    % set(gcf,'Name',['Global solution H_U over fictitious domain'])
elseif ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Global solution U_' num2str(i) ' over fictitious domain'])
    % set(gcf,'Name',['Global solution U_' num2str(i) ' over fictitious domain'])
elseif ischarin('rotation',varargin)
    i = getcharin('rotation',varargin);
    figure('Name',['Global solution rot_U_' num2str(i) ' over fictitious domain'])
    % set(gcf,'Name',['Global solution rot_U_' num2str(i) ' over fictitious domain'])
else
    figure('Name','Global solution U over fictitious domain')
    % set(gcf,'Name','Global solution U over fictitious domain')
end
clf
set(gcf,'Color','w')

Ut = setevolparam(Ut,'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,...
    'view',p.Results.view,'camup',p.Results.camup,'campos',p.Results.campos,'FontSize',p.Results.FontSize);
frame = evol_sol(Ut,glob.S,'rescale',p.Results.rescale,varargin{:}); % save the frames

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
