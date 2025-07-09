function mov = evolSolutionCell(T,St,ut,varargin)
% function mov = evolSolutionCell(T,St,ut,varargin)
% Display evolution of solution ut
% T : TIMEMODEL
% St: MODEL or cell array of MODEL
% ut: cell array of solution u
% mov: movie

p = ImprovedInputParser;
addParameter(p,'rescale',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'view',[],@isnumeric);
addParameter(p,'camup','auto',@(x) isnumeric(x) || ischar(x));
addParameter(p,'campos','auto',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
addParameter(p,'axison',false,@islogical);
addParameter(p,'boxon',false,@islogical);
addParameter(p,'boxstylefull',false,@islogical);
addParameter(p,'noxtick',false,@islogical);
addParameter(p,'noytick',false,@islogical);
addParameter(p,'noztick',false,@islogical);
addParameter(p,'plotiter',false,@islogical);
addParameter(p,'plottime',true,@islogical);
addParameter(p,'filename','solution',@ischar);
addParameter(p,'pathname','./',@ischar);
addParameter(p,'formats',{'avi','mp4'},@(x) ischar(x) || iscell(x));
addParameter(p,'FrameRate',30,@isnumeric);
addParameter(p,'Quality',100,@isnumeric);
parse(p,varargin{:})

varargin = delcharin({'rescale','colorbar','colormap','view','camup','campos','FontSize',...
    'axison','boxon','boxstylefull','noxtick','noytick','noztick',...
    'plotiter''plottime','filename','pathname','formats','FrameRate','Quality'},varargin);
if isa(p.Results.formats,'char')
    p.Results.formats = {p.Results.formats};
end

% Create a new figure
if ischarin('sigma',varargin)
    i = getcharin('sigma',varargin);
    figure('Name',['Solution sig_' num2str(i)])
    % set(gcf,'Name',['Solution sig_' num2str(i)])
elseif ischarin('epsilon',varargin)
    i = getcharin('epsilon',varargin);
    figure('Name',['Solution eps_' num2str(i)])
    % set(gcf,'Name',['Solution eps_' num2str(i)])
elseif ischarin('energyint',varargin)
    figure('Name',['Solution H'])
    % set(gcf,'Name',['Solution H'])
elseif ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Solution u_' num2str(i)])
    % set(gcf,'Name',['Solution u_' num2str(i)])
elseif ischarin('rotation',varargin)
    i = getcharin('rotation',varargin);
    figure('Name',['Solution r_' num2str(i)])
    % set(gcf,'Name',['Solution r_' num2str(i)])
else
    figure('Name','Solution u')
    % set(gcf,'Name','Solution u')
end
clf
set(gcf,'color','w')

T = setevolparam(T,'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,...
    'view',p.Results.view,'camup',p.Results.camup,'campos',p.Results.campos,'FontSize',p.Results.FontSize,...
    'axison',p.Results.axison,'boxon',p.Results.boxon,'boxstylefull',p.Results.boxstylefull,...
    'noxtick',p.Results.noxtick,'noytick',p.Results.noytick,'noztick',p.Results.noztick,...
    'plotiter',p.Results.plotiter,'plottime',p.Results.plottime);
frame = evol_sol_cell(T,ut,St,'rescale',p.Results.rescale,varargin{:}); % save the frames

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
