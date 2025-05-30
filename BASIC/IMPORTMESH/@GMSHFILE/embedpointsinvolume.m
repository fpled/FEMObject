function u = embedpointsinvolume(u,numberpoints,numbervolume)
% function u = embedpointsinvolume(u,numberpoints,numbervolume)

for k=1:length(numberpoints)
    u = embedpointinvolume(u,numberpoints(k),numbervolume);
end
