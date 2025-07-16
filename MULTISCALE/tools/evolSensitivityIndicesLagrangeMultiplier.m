function mov = evolSensitivityIndicesLagrangeMultiplier(interfaces,T,lambdat,alpha,varargin)
% function mov = evolSensitivityIndicesLagrangeMultiplier(interfaces,T,lambdat,alpha,varargin)
% Display evolution of the sensitivity indices of Lagrange multiplier lambdat
% associated with the group of variables alpha in {1,..,d} with d = ndims(lambdat)
% interfaces: Interfaces or Interface
% T: TIMEMODEL
% lambdat: FunctionalBasis of Lagrange multiplier lambda
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

varargin = delcharin({'rescale','colorbar','colormap','view','camup','campos','FontSize',....
    'filename','pathname','formats','FrameRate','Quality'},varargin);
if isa(p.Results.formats,'char')
    p.Results.formats = {p.Results.formats};
end

if isa(interfaces,'Interfaces')
    numbers = getnumber(interfaces);
    if ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Sensitivity index of lambda_' num2str(i) ' over interfaces #' num2str([numbers{:}]) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Sensitivity index of lambda_' num2str(i) ' over interfaces #' num2str([numbers{:}]) ' for random variables #' num2str(alpha)])
    else
        figure('Name',['Sensitivity index of lambda over interfaces #' num2str([numbers{:}]) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Sensitivity index of lambda over interfaces #' num2str([numbers{:}]) ' for random variables #' num2str(alpha)])
    end
    clf
    set(gcf,'Color','w')
    
    interface = interfaces.interfaces;
    S = cellfun(@(interface) interface.S,interface,'UniformOutput',false);
    slambdat = cellfun(@(interface) reshape(varianceConditionalExpectation(lambdat{interface.number},alpha)./max(variance(lambdat{interface.number})),lambdat{interface.number}.sz),interface,'UniformOutput',false);
    slambdat = cellfun(@(s) TIMEMATRIX(s,T),slambdat,'UniformOutput',false);
    frame = evol_sol(T,slambdat,S,'rescale',p.Results.rescale,varargin{:}); % save the frames
elseif isa(interfaces,'Interface')
    interface = interfaces;
    if ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Sensitivity index of Lagrange multiplier lambda_' num2str(i) ' over interface #' num2str(interface.number) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Sensitivity index of Lagrange multiplier lambda_' num2str(i) ' over interface #' num2str(interface.number) ' for random variables #' num2str(alpha)])
    else
        figure('Name',['Sensitivity index of Lagrange multiplier lambda over interface #' num2str(interface.number) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Sensitivity index of Lagrange multiplier lambda over interface #' num2str(interface.number) ' for random variables #' num2str(alpha)])
    end
    clf
    set(gcf,'Color','w')
    
    sz = lambdat{interface.number}.sz;
    vlambdat = variance(lambdat{interface.number});
    slambdat = varianceConditionalExpectation(lambdat{interface.number},alpha)./max(vlambdat);
    slambdat = reshape(slambdat,sz);
    slambdat = TIMEMATRIX(slambdat,T);
    
    slambdat = setevolparam(slambdat,'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,...
        'view',p.Results.view,'camup',p.Results.camup,'campos',p.Results.campos,'FontSize',p.Results.FontSize);
    frame = evol_sol(slambdat,interface.S,'rescale',p.Results.rescale,varargin{:}); % save the frames
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
