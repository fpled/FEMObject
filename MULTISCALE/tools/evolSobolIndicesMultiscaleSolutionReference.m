function mov = evolSobolIndicesMultiscaleSolutionReference(glob,patches,interfaces,T,Ut_ref,wt_ref,alpha,varargin)
% function mov = evolSobolIndicesMultiscaleSolutionReference(glob,patches,interfaces,T,Ut_ref,wt_ref,alpha,varargin)
% Display evolution of the Sobol indices of reference multiscale solution ut_ref=(Ut_ref,wt_ref)
% associated with the group of variables alpha in {1,..,d} with d = ndims(ut_ref)
% glob: Global or GlobalOutside
% patches: Patches
% interfaces: Interfaces
% T: TIMEMODEL
% Ut_ref: FunctionalBasis of reference global solution U
% wt_ref: FunctionalBasis of reference local solution w
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
n = numel(patches);

if ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Sobol index of u_ref_' num2str(i) '=(U_ref_' num2str(i) ',w_ref_' num2str(i) ') over domain for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Sobol index of u_ref_' num2str(i) '=(U_ref_' num2str(i) ',w_ref_' num2str(i) ') over domain for random variables #' num2str(alpha)])
else
    figure('Name',['Sobol index of u_ref=(U_ref,w_ref) over domain for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Sobol index of u_ref=(U_ref,w_ref) over domain for random variables #' num2str(alpha)])
end
clf
set(gcf,'Color','w')

if isa(glob,'Global')
    S_out = glob.S_out;
elseif isa(glob,'GlobalOutside')
    S_out = glob.S;
end

sz = Ut_ref.sz;
dUt_ref = ndims(Ut_ref);
sUt_ref = SensitivityAnalysis.sobolIndices(Ut_ref,alpha,dU_ref);
sUt_ref = reshape(sUt_ref,sz);
sUt_ref = TIMEMATRIX(sUt_ref,T);
sUt_ref = unfreevector(S_out,sUt_ref)-calc_init_dirichlet(S_out)*one(T);

patch = patches.patches;
interface = interfaces.interfaces;
S_patch = cellfun(@(patch) patch.S,patch,'UniformOutput',false);
S_interface = cellfun(@(interface) interface.S,interface,'UniformOutput',false);
swt_ref = cellfun(@(patch) TIMEMATRIX(reshape(SensitivityAnalysis.sobolIndices(wt_ref{patch.number},alpha,ndims(wt_ref{patch.number})),wt_ref{patch.number}.sz),T),patch,'UniformOutput',false);
slambdat_ref = cellfun(@(patch,interface) interface.P_patch*swt_ref{patch.number},patch,interface,'UniformOutput',false);

T = setevolparam(T,'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,...
    'view',p.Results.view,'camup',p.Results.camup,'campos',p.Results.campos,'FontSize',p.Results.FontSize);
frame = evol_sol(T,[{sUt_ref},swt_ref,slambdat_ref],[{S_out},S_patch,S_interface],'rescale',p.Results.rescale,'interfaces',2+n:1+2*n,varargin{:}); % save the frames

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
