function mov = evolMeanLagrangeMultiplier(interfaces,T,lambdat,varargin)
% function mov = evolMeanLagrangeMultiplier(interfaces,T,lambdat,varargin)
% Display evolution of the mean (mathematical expectation) of Lagrange multiplier lambdat
% interfaces: Interfaces or Interface
% T: TIMEMODEL
% lambdat: FunctionalBasis of Lagrange multiplier lambda
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

if isa(interfaces,'Interfaces')
    numbers = getnumber(interfaces);
    if ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Mean of lambda_' num2str(i) ' over interfaces #' num2str([numbers{:}])])
        % set(gcf,'Name',['Mean of lambda_' num2str(i) ' over interfaces #' num2str([numbers{:}])])
    else
        figure('Name',['Mean of lambda over interfaces #' num2str([numbers{:}])])
        % set(gcf,'Name',['Mean of lambda over interfaces #' num2str([numbers{:}])])
    end
    clf
    set(gcf,'Color','w')
    
    interface = interfaces.interfaces;
    S = cellfun(@(interface) interface.S,interface,'UniformOutput',false);
    lambdat = cellfun(@(interface) TIMEMATRIX(reshape(mean(lambdat{interface.number}),lambdat{interface.number}.sz),T),interface,'UniformOutput',false);
    frame = evol_sol(T,lambdat,S,'rescale',p.Results.rescale,varargin{:}); % save the frames
elseif isa(interfaces,'Interface')
    interface = interfaces;
    if ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Mean of Lagrange multiplier lambda_' num2str(i) ' over interface #' num2str(interface.number)])
        % set(gcf,'Name',['Mean of Lagrange multiplier lambda_' num2str(i) ' over interface #' num2str(interface.number)])
    else
        figure('Name',['Mean of Lagrange multiplier lambda over interface #' num2str(interface.number)])
        % set(gcf,'Name',['Mean of Lagrange multiplier lambda over interface #' num2str(interface.number)])
    end
    clf
    set(gcf,'Color','w')
    
    sz = lambdat{interface.number}.sz;
    lambdat{interface.number} = mean(lambdat{interface.number});
    lambdat{interface.number} = reshape(lambdat{interface.number},sz);
    lambdat{interface.number} = TIMEMATRIX(lambdat{interface.number},T);
    
    lambdat{interface.number} = setevolparam(lambdat{interface.number},'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,...
        'view',p.Results.view,'camup',p.Results.camup,'campos',p.Results.campos,'FontSize',p.Results.FontSize);
    frame = evol_sol(lambdat{interface.number},interface.S,'rescale',p.Results.rescale,varargin{:}); % save the frames
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
