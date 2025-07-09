function mov = evolGlobalSolutionReference(glob,Ut_ref,varargin)
% function mov = evolGlobalSolutionReference(glob,Ut_ref,varargin)
% Display evolution of reference global solution Ut_ref
% glob: Global or GlobalOutside
% Ut_ref: TIMEMATRIX of reference global solution U
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
    figure('Name',['Reference global solution sig_U_ref_' num2str(i) ' over complementary subdomain'])
    % set(gcf,'Name',['Reference global solution sig_U_ref_' num2str(i) ' over complementary subdomain'])
elseif ischarin('epsilon',varargin)
    i = getcharin('epsilon',varargin);
    figure('Name',['Reference global solution eps_U_ref_' num2str(i) ' over complementary subdomain'])
    % set(gcf,'Name',['Reference global solution eps_U_ref_' num2str(i) ' over complementary subdomain'])
elseif ischarin('energyint',varargin)
    figure('Name',['Reference global solution H_U_ref over complementary subdomain'])
    % set(gcf,'Name',['Reference global solution H_U_ref over complementary subdomain'])
elseif ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Reference global solution U_ref_' num2str(i) ' over complementary subdomain'])
    % set(gcf,'Name',['Reference global solution U_ref_' num2str(i) ' over complementary subdomain'])
elseif ischarin('rotation',varargin)
    i = getcharin('rotation',varargin);
    figure('Name',['Reference global solution rot_U_ref_' num2str(i) ' over complementary subdomain'])
    % set(gcf,'Name',['Reference global solution rot_U_ref_' num2str(i) ' over complementary subdomain'])
else
    figure('Name','Reference global solution U_ref over complementary subdomain')
    % set(gcf,'Name','Reference global solution U_ref over complementary subdomain')
end
clf
set(gcf,'color','w')

if isa(glob,'Global')
    S_out = glob.S_out;
elseif isa(glob,'GlobalOutside')
    S_out = glob.S;
end

Ut_ref = setevolparam(Ut_ref,'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,...
    'view',p.Results.view,'camup',p.Results.camup,'campos',p.Results.campos,'FontSize',p.Results.FontSize);
frame = evol_sol(Ut_ref,S_out,'rescale',p.Results.rescale,varargin{:}); % save the frames

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
