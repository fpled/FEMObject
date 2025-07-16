function varargout = plotenveloppe(Y,face,varargin)
% function plotenveloppe(Y,facecolor)
% Y : TIMEMATRIX
% facecolor : couleur de remplissage
% 
% function plotenveloppe(Y,facecolor,'EdgeColor',edgecolor)
% edgecolor : couleurs des contrours ('k' par defaut)
%
% function plotenveloppe(Y,facecolor,'FaceAlpha',alpha,'EdgeAlpha',edgealpha)
% alpha et edgealpha : coefficients de transparence (0.5 par defaut)

if nargin<=1
 face = 'y';
end

[t,rep]=gettplot(Y);

if ~isa(Y.value,'double')
    error('pas programme')
end

if all(size(Y)==1)
Ymax = max(Y,0);
Ymin = -max(-Y,0);
else
Ymin = min(Y,[],1);
Ymax = max(Y,[],1);
end

edge = getcharin('EdgeColor',varargin,'k');
transp = getcharin('FaceAlpha',varargin,0.5);
transpedge = getcharin('EdgeAlpha',varargin,transp);

t=[t,fliplr(t)];
Y = [Ymin.value(:,rep),fliplr(Ymax.value(:,rep))];
fillhandle = fill(t,Y,face);

set(fillhandle,'EdgeColor',edge,'FaceAlpha',transp,'EdgeAlpha',transpedge);

if nargout==1
  varargout{1}=fillhandle;  
end