function mov = evolModel(T,St,varargin)
% function mov = evolSolution(T,St,varargin)
% Display evolution of model St
% T : TIMEMODEL
% St: cell array of MODEL
% mov: movie

p = ImprovedInputParser;
addParameter(p,'legend',true,@islogical);
addParameter(p,'node',false,@isscalar);
addParameter(p,'FontSize',16,@isscalar);
addParameter(p,'LineWidth',0.5,@isscalar);
addParameter(p,'Interpreter','latex',@ischar);
addParameter(p,'FaceColor',getfacecolor(1),@(x) isnumeric(x) || ischar(x));
parse(p,varargin{:})

varargin = delcharin({'legend','node','FontSize','LineWidth','Interpreter','FaceColor'},varargin);

p = ImprovedInputParser;
addParameter(p,'legend',true,@islogical);
addParameter(p,'node',false,@(x) islogical(x) || ischar(x));
addParameter(p,'view',[],@isnumeric);
addParameter(p,'camup','auto',@(x) isnumeric(x) || ischar(x));
addParameter(p,'campos','auto',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
addParameter(p,'LineWidth',0.5,@isscalar);
addParameter(p,'Interpreter','latex',@ischar);
addParameter(p,'FaceColor',getfacecolor(1),@(x) isnumeric(x) || ischar(x));
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

varargin = delcharin({'view','camup','campos','FontSize','LineWidth','Interpreter','FaceColor',...
    'axison','boxon','boxstylefull','noxtick','noytick','noztick',...
    'plotiter''plottime','filename','pathname','formats','FrameRate','Quality'},varargin);
if isa(p.Results.formats,'char')
    p.Results.formats = {p.Results.formats};
end

figure('Name','Mesh')
% set(gcf,'Name','Mesh')
clf
set(gcf,'Color','w')

T = setevolparam(T,'view',p.Results.view,'camup',p.Results.camup,'campos',p.Results.campos,'FontSize',p.Results.FontSize,...
    'axison',p.Results.axison,'boxon',p.Results.boxon,'boxstylefull',p.Results.boxstylefull,...
    'noxtick',p.Results.noxtick,'noytick',p.Results.noytick,'noztick',p.Results.noztick,'plotiter',p.Results.plotiter,'plottime',p.Results.plottime);
frame = evol_model(T,St,varargin{:}); % save the frames

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
