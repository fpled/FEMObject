function mov = evolMultiscaleSolution(glob,patches,interfaces,Ut,wt,varargin)
% function mov = evolMultiscaleSolution(glob,patches,interfaces,Ut,wt,varargin)
% Display evolution of multiscale solution ut=(Ut,wt)
% glob: Global
% patches: Patches
% interfaces: Interfaces
% Ut: TIMEMATRIX of global solution U
% wt: TIMEMATRIX of local solution w
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
    figure('Name',['Multiscale solution sig_u_' num2str(i) '=(sig_U_' num2str(i) ',sig_w_' num2str(i) ') over domain'])
    % set(gcf,'Name',['Multiscale solution sig_u_' num2str(i) '=(sig_U_' num2str(i) ',sig_w_' num2str(i) ') over domain'])
elseif ischarin('epsilon',varargin)
    i = getcharin('epsilon',varargin);
    figure('Name',['Multiscale solution eps_u_' num2str(i) '=(eps_U_' num2str(i) ',eps_w_' num2str(i) ') over domain'])
    % set(gcf,'Name',['Multiscale solution eps_u_' num2str(i) '=(eps_U_' num2str(i) ',eps_w_' num2str(i) ') over domain'])
elseif ischarin('energyint',varargin)
    figure('Name',['Multiscale solution H_u=(H_U,H_w) over domain'])
    % set(gcf,'Name',['Multiscale solution H_u=(H_U,H_w) over domain'])
elseif ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Multiscale solution u_' num2str(i) '=(U_' num2str(i) ',w_' num2str(i) ') over domain'])
    % set(gcf,'Name',['Multiscale solution u_' num2str(i) '=(U_' num2str(i) ',w_' num2str(i) ') over domain'])
elseif ischarin('rotation',varargin)
    i = getcharin('rotation',varargin);
    figure('Name',['Multiscale solution rot_u_' num2str(i) '=(rot_U_' num2str(i) ',rot_w_' num2str(i) ') over domain'])
    % set(gcf,'Name',['Multiscale solution rot_u_' num2str(i) '=(rot_U_' num2str(i) ',rot_w_' num2str(i) ') over domain'])
else
    figure('Name','Multiscale solution u=(U,w) over domain')
    % set(gcf,'Name','Multiscale solution u=(U,w) over domain')
end
clf
set(gcf,'color','w')

if size(Ut,1)==glob.S.nbddl
    P_out = calcProjection(glob,'free',false);
else
    P_out = glob.P_out;
end

patch = patches.patches;
interface = interfaces.interfaces;
S_patch = cellfun(@(patch) patch.S,patch,'UniformOutput',false);
S_interface = cellfun(@(interface) interface.S,interface,'UniformOutput',false);
lambdat = cellfun(@(patch,interface) interface.P_patch*wt{patch.number},patch,interface,'UniformOutput',false);

T = gettimemodel(Ut);
T = setevolparam(T,'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,...
    'view',p.Results.view,'camup',p.Results.camup,'campos',p.Results.campos,'FontSize',p.Results.FontSize);
frame = evol_sol(T,[{P_out*Ut},wt,lambdat],[{glob.S_out},S_patch,S_interface],'rescale',p.Results.rescale,'interfaces',2+n:1+2*n,varargin{:}); % save the frames

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
