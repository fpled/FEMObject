function mov = evolClosedSobolIndicesLagrangeMultiplierReference(interfaces,T,lambdat_ref,alpha,varargin)
% function mov = evolClosedSobolIndicesLagrangeMultiplierReference(interfaces,T,lambdat_ref,alpha,varargin)
% Display evolution of the Closed Sobol indices of reference Lagrange multiplier lambdat_ref
% associated with the group of variables alpha in {1,..,d} with d = ndims(lambdat)
% interfaces: Interfaces or Interface
% T: TIMEMODEL
% lambdat_ref: FunctionalBasis of reference Lagrange multiplier lambda
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

if isa(interfaces,'Interfaces')
    numbers = getnumber(interfaces);
    if ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Sobol index of lambda_ref_' num2str(i) ' over interfaces #' num2str([numbers{:}]) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Sobol index of lambda_ref_' num2str(i) ' over interfaces #' num2str([numbers{:}]) ' for random variables #' num2str(alpha)])
    else
        figure('Name',['Sobol index of lambda_ref over interfaces #' num2str([numbers{:}]) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Sobol index of lambda_ref over interfaces #' num2str([numbers{:}]) ' for random variables #' num2str(alpha)])
    end
    clf
    set(gcf,'Color','w')
    
    interface = interfaces.interfaces;
    S = cellfun(@(interface) interface.S,interface,'UniformOutput',false);
    slambdat_ref = cellfun(@(interface) TIMEMATRIX(reshape(SensitivityAnalysis.closedSobolIndices(lambdat_ref{interface.number},alpha,ndims(lambdat_ref{interface.number})),lambdat_ref{interface.number}.sz),T),interface,'UniformOutput',false);
    frame = evol_sol(T,slambdat_ref,S,'rescale',p.Results.rescale,varargin{:}); % save the frames
elseif isa(interfaces,'Interface')
    interface = interfaces;
    if ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Closed Sobol index of Lagrange multiplier lambda_ref_' num2str(i) ' over interface #' num2str(interface.number) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Closed Sobol index of Lagrange multiplier lambda_ref_' num2str(i) ' over interface #' num2str(interface.number) ' for random variables #' num2str(alpha)])
    else
        figure('Name',['Closed Sobol index of Lagrange multiplier lambda_ref over interface #' num2str(interface.number) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Closed Sobol index of Lagrange multiplier lambda_ref over interface #' num2str(interface.number) ' for random variables #' num2str(alpha)])
    end
    clf
    set(gcf,'Color','w')
    
    sz = lambdat_ref{interface.number}.sz;
    slambdat_ref = SensitivityAnalysis.closedSobolIndices(lambdat_ref{interface.number},alpha);
    slambdat_ref = reshape(slambdat_ref,sz);
    slambdat_ref = TIMEMATRIX(slambdat_ref,T);
    
    slambdat_ref = setevolparam(slambdat_ref,'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,...
        'view',p.Results.view,'camup',p.Results.camup,'campos',p.Results.campos,'FontSize',p.Results.FontSize);
    frame = evol_sol(slambdat_ref,interface.S,'rescale',p.Results.rescale,varargin{:}); % save the frames
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
