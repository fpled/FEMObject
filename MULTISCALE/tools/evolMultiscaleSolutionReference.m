function mov = evolMultiscaleSolutionReference(glob,patches,interfaces,Ut_ref,wt_ref,varargin)
% function mov = evolMultiscaleSolutionReference(glob,patches,interfaces,Ut_ref,wt_ref,varargin)
% Display evolution of reference multiscale solution ut_ref=(Ut_ref,wt_ref)
% glob: Global or GlobalOutside
% patches: Patches
% interfaces: Interfaces
% Ut_ref: TIMEMATRIX of reference global solution U
% wt_ref: TIMEMATRIX of reference local solution w
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
n = numel(patches);

if ischarin('sigma',varargin)
    i = getcharin('sigma',varargin);
    figure('Name',['Reference multiscale solution sig_u_ref_' num2str(i) '=(sig_U_ref_' num2str(i) ',sig_w_ref_' num2str(i) ') over domain'])
    % set(gcf,'Name',['Reference multiscale solution sig_u_ref_' num2str(i) '=(sig_U_ref_' num2str(i) ',sig_w_ref_' num2str(i) ') over domain'])
elseif ischarin('epsilon',varargin)
    i = getcharin('epsilon',varargin);
    figure('Name',['Reference multiscale solution eps_u_ref_' num2str(i) '=(eps_U_ref_' num2str(i) ',eps_w_ref_' num2str(i) ') over domain'])
    % set(gcf,'Name',['Reference multiscale solution eps_u_ref_' num2str(i) '=(eps_U_ref_' num2str(i) ',eps_w_ref_' num2str(i) ') over domain'])
elseif ischarin('energyint',varargin)
    figure('Name',['Reference multiscale solution H_u_ref=(H_U_ref,H_w_ref) over domain'])
    % set(gcf,'Name',['Reference multiscale solution H_u_ref=(H_U_ref,H_w_ref) over domain'])
elseif ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Reference multiscale solution u_ref_' num2str(i) '=(U_ref_' num2str(i) ',w_ref_' num2str(i) ') over domain'])
    % set(gcf,'Name',['Reference multiscale solution u_ref_' num2str(i) '=(U_ref_' num2str(i) ',w_ref_' num2str(i) ') over domain'])
elseif ischarin('rotation',varargin)
    i = getcharin('rotation',varargin);
    figure('Name',['Reference multiscale solution rot_u_ref_' num2str(i) '=(rot_U_ref_' num2str(i) ',rot_w_ref_' num2str(i) ') over domain'])
    % set(gcf,'Name',['Reference multiscale solution rot_u_ref_' num2str(i) '=(rot_U_ref_' num2str(i) ',rot_w_ref_' num2str(i) ') over domain'])
else
    figure('Name','Reference multiscale solution u_ref=(U_ref,w_ref) over domain')
    % set(gcf,'Name','Reference multiscale solution u_ref=(U_ref,w_ref) over domain')
end
clf
set(gcf,'Color','w')

if isa(glob,'Global')
    S_out = glob.S_out;
elseif isa(glob,'GlobalOutside')
    S_out = glob.S;
end

patch = patches.patches;
interface = interfaces.interfaces;
S_patch = cellfun(@(patch) patch.S,patch,'UniformOutput',false);
S_interface = cellfun(@(interface) interface.S,interface,'UniformOutput',false);
lambdat_ref = cellfun(@(patch,interface) interface.P_patch*wt_ref{patch.number},patch,interface,'UniformOutput',false);

T = gettimemodel(Ut_ref);
T = setevolparam(T,'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,...
    'view',p.Results.view,'camup',p.Results.camup,'campos',p.Results.campos,'FontSize',p.Results.FontSize);
frame = evol_sol(T,[{Ut_ref},wt_ref,lambdat_ref],[{S_out},S_patch,S_interface],'rescale',p.Results.rescale,'interfaces',2+n:1+2*n,varargin{:}); % save the frames

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
