function u = embedcurvesinvolume(u,numberlines,numbervolume)
% function u = embedcurvesinvolume(u,numberlines,numbervolume)

for k=1:length(numberlines)
    u = embedcurveinvolume(u,numberlines(k),numbervolume);
end
