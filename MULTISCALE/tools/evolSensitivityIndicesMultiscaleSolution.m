function mov = evolSensitivityIndicesMultiscaleSolution(glob,patches,interfaces,T,Ut,wt,alpha,varargin)
% function mov = evolSensitivityIndicesMultiscaleSolution(glob,patches,interfaces,T,Ut,wt,alpha,varargin)
% Display evolution of the sensitivity indices of multiscale solution ut=(Ut,wt)
% associated with the group of variables alpha in {1,..,d} with d = ndims(ut)
% glob: Global
% patches: Patches
% interfaces: Interfaces
% T: TIMEMODEL
% Ut: FunctionalBasis of global solution U
% wt: FunctionalBasis of local solution w
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
    figure('Name',['Sensitivity index of u_' num2str(i) '=(U_' num2str(i) ',w_' num2str(i) ') over domain for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Sensitivity index of u_' num2str(i) '=(U_' num2str(i) ',w_' num2str(i) ') over domain for random variables #' num2str(alpha)])
else
    figure('Name',['Sensitivity index of u=(U,w) over domain for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Sensitivity index of u=(U,w) over domain for random variables #' num2str(alpha)])
end
clf
set(gcf,'color','w')

sz = Ut.sz;
vUt = variance(Ut);
vUt = reshape(vUt,sz);
vUt = TIMEMATRIX(vUt,T);

if size(vUt,1)==glob.S.nbddl
    P_out = calcProjection(glob,'free',false);
else
    P_out = glob.P_out;
end
vUt_out = P_out*vUt;

patch = patches.patches;
interface = interfaces.interfaces;
vwt = cellfun(@(patch) TIMEMATRIX(reshape(variance(wt{patch.number}),wt{patch.number}.sz),T),patch,'UniformOutput',false);
max_var = getvalue(max(vUt_out));
for k=1:n
    max_var = max(max_var,getvalue(max(vwt{k})));
end

vUt = varianceConditionalExpectation(Ut,alpha);
vUt = reshape(vUt,sz);
vUt_out = P_out*vUt;
sUt_out = vUt_out./repmat(max_var,size(vUt_out,1),1);
sUt_out = TIMEMATRIX(sUt_out,T);
sUt_out = unfreevector(glob.S_out,sUt_out)-calc_init_dirichlet(glob.S_out)*one(T);

S_patch = cellfun(@(patch) patch.S,patch,'UniformOutput',false);
S_interface = cellfun(@(interface) interface.S,interface,'UniformOutput',false);
vwt = cellfun(@(patch) reshape(varianceConditionalExpectation(wt{patch.number},alpha),wt{patch.number}.sz),patch,'UniformOutput',false);
swt = cellfun(@(vwt) vwt./repmat(max_var,size(vwt,1),1),vwt,'UniformOutput',false);
swt = cellfun(@(swt) TIMEMATRIX(swt,T),swt,'UniformOutput',false);
slambdat = cellfun(@(patch,interface) interface.P_patch*swt{patch.number},patch,interface,'UniformOutput',false);

T = setevolparam(T,'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,...
    'view',p.Results.view,'camup',p.Results.camup,'campos',p.Results.campos,'FontSize',p.Results.FontSize);
frame = evol_sol(T,[{sUt_out},swt,slambdat],[{glob.S_out},S_patch,S_interface],'rescale',p.Results.rescale,'interfaces',2+n:1+2*n,varargin{:}); % save the frames

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
