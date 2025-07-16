function mov = evolLagrangeMultiplier(interfaces,lambdat,varargin)
% function mov = evolLagrangeMultiplier(interfaces,lambdat,varargin)
% Display evlution of Lagrange multiplier lambda
% interfaces: Interfaces or Interface
% lambdat: TIMEMATRIX of Lagrange multiplier lambda
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
        figure('Name',['Lagrange mulitplier lambda_' num2str(i) ' over interfaces #' num2str([numbers{:}])])
        % set(gcf,'Name',['Lagrange mulitplier lambda_' num2str(i) ' over interfaces #' num2str([numbers{:}])])
    else
        figure('Name',['Lagrange mulitplier lambda over interfaces #' num2str([numbers{:}])])
        % set(gcf,'Name',['Lagrange mulitplier lambda over interfaces #' num2str([numbers{:}])])
    end
    clf
    set(gcf,'Color','w')
    
    n = numel(interfaces);
    S = cellfun(@(interface) interface.S,interfaces.interfaces,'UniformOutput',false);
    for k=1:n
        interface = interfaces.interfaces{k};
        lambdat{interface.number} = setevolparam(lambdat{interface.number},'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,...
            'view',p.Results.view,'camup',p.Results.camup,'campos',p.Results.campos,'FontSize',p.Results.FontSize);
    end
    T = gettimemodel(lambdat{1});
    frame = evol_sol(T,lambdat,S,'rescale',p.Results.rescale,varargin{:}); % save the frames
elseif isa(interfaces,'Interface')
    interface = interfaces;
    if ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Lagrange mulitplier lambda_' num2str(i) ' over interface #' num2str(interface.number)])
        % set(gcf,'Name',['Lagrange mulitplier lambda_' num2str(i) ' over interface #' num2str(interface.number)])
    else
        figure('Name',['Lagrange mulitplier lambda over interface #' num2str(interface.number)])
        % set(gcf,'Name',['Lagrange mulitplier lambda over interface #' num2str(interface.number)])
    end
    clf
    set(gcf,'Color','w')
    
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

