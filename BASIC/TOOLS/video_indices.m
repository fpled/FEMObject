function video_indices(PC,varargin)
% function video_indices(PC,'dim',dim,'filename',filename,'pathname',pathname)
% Make a video of the adaptive construction of the multi-index set of
% POLYCHAOS PC projected on dimension dim using adaptive algorithm
% PC: cell of POLYCHAOS
% dim: dimension, 1:min(3,getM(PC{end})) by default
% filename: movie filename, 'mutli_index_set.avi' by default
% pathname: movie pathname, './' by default

p = inputParser;
addParameter(p,'dim',1:min(3,getM(PC{end})),@isnumeric);
addParameter(p,'pause',false,@islogical);
addParameter(p,'filename','multi_index_set',@ischar);
addParameter(p,'pathname','./',@ischar);
parse(p,varargin{:})

[ok,rep]=ismember(p.Results.dim,PC{end});
if ~all(ok)
    error('Wrong dimensions')
end

% Create a new figure and save position
h = figure('Name','Multi-index set');
clf
set(gcf,'color','w')
winsize = get(h,'Position');
winsize(1:2) = [0 0]; % set borders of the movie
% Create movie file
filename = [p.Results.filename '_dim' sprintf('_%d',p.Results.dim(1:end))];
mov = VideoWriter(fullfile(p.Results.pathname,filename),'MPEG-4');
mov.FrameRate = 30;
mov.Quality = 100;
open(mov);
% set(h,'NextPlot','replacechildren'); % ensure each frame is of the same size as before, so all frames will appear in the same movie

PC_dim = projectdim(PC{end},p.Results.dim);
ind = getindices(PC_dim);

switch length(p.Results.dim)
    case 1
        ax = [-eps max(ind(:,1))+1 -eps 1];
        tick = 0:max(ind(:,1));
    case 2
        ax = [-eps max(ind(:,1))+1 -eps max(ind(:,2))+1];
        tick = [{0:max(ind(:,1))} {0:max(ind(:,2))}];
    case 3
        ax = [-eps max(ind(:,1))+1 -eps max(ind(:,2))+1 -eps max(ind(:,3))+1];
        tick = [{0:max(ind(:,1))} {0:max(ind(:,2))} {0:max(ind(:,3))}];
end
for i=1:numel(PC)
    if i==1
        plot_indices(PC{i},'dim',p.Results.dim,'legend',false,'axis',ax,'tick',tick)
        hold on
        frame = getframe(h,winsize); % save the frame
        writeVideo(mov,frame); % add this frame to the movie
        if p.Results.pause
            pause
        end
    else
        if ~isequal(PC{i},PC{i-1})
            PC_add = setdiff(PC{i},PC{i-1});
            plot_indices(PC_add,'dim',p.Results.dim,'MarkerColor','r','legend',false,'axis',ax,'tick',tick)
            hold off
            frame = getframe(h,winsize); % save the frame
            writeVideo(mov,frame); % add this frame to the movie
            if p.Results.pause
                pause
            end
            plot_indices(PC{i},'dim',p.Results.dim,'MarkerColor','b','legend',false,'axis',ax,'tick',tick)
            hold on
            frame = getframe(h,winsize); % save the frame
            writeVideo(mov,frame); % add this frame to the movie
        end
    end
end
hold off
close(mov);

end
