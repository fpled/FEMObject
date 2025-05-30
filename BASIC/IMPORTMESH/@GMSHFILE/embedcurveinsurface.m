function u = embedcurveinsurface(u,numberline,numbersurface)
% function u = embedcurveinsurface(u,numberline,numbersurface)

u = embed(u,'Curve',numberline,'Surface',numbersurface);
