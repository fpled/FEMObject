function u = embedcurvesinsurface(u,numberlines,numbersurface)
% function u = embedcurvesinsurface(u,numberlines,numbersurface)

for k=1:length(numberlines)
    u = embedcurveinsurface(u,numberlines(k),numbersurface);
end
