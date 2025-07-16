function varargout = sampleplot(us,varargin)
% function [P,x,y] = samleplot(us,varargin)
% trace la pdf ou cdf d'un ensemble de realisations
% plot(us,courbe,varargin) avec courbe  = 'pdf' ou 'cdf' ou 'multipdf'(pdf par defaut)
% pdfplot(us,'npts',n) trace la courbe avec n points (100 par defaut)
% pdfplot(us,'bar') trace des barres.
% pdfplot(us,'nbbar',b) facteur multiplicatif pour le nombre de points (1 par defaut)
% tous les autres arguments sont passes en argument de plot ou bar

if any(size(us)==1)
    us=us(:);
end
nbbar=getcharin('nbbar',varargin,1);
N = size(us,1);
n=getcharin('npts',varargin,floor(1+3.3*log(N)/log(10))*3*nbbar+1);
tracebar=ischarin('bar',varargin);
varargin=delonlycharin('bar',varargin);
varargin = delcharin({'npts','nbs','nbbar'},varargin);
[rep,post]=ischarin({'pdf','cdf','icdf','multipdf'},varargin);

withks =ischarin('ksdensity',varargin);
varargin=delonlycharin('ksdensity',varargin);

if ~any(rep)
    courbe = 'pdf'    ;
else
    courbe = varargin{post(find(rep))};
end
varargin = delonlycharin({'pdf','cdf','icdf','multipdf'},varargin);

switch courbe
case 'cdf'
    [P,x] = ecdf(us);
case 'pdf'
    if ischarin('axis',varargin)
        ax = getcharin('axis',varargin);
        varargin = delcharin('axis',varargin);
        x=linspace(ax(1),ax(2),n);
    else
        x=linspace(min(us)-10*eps,max(us)+10*eps,n);
    end
    if withks
        [P,x] = ksdensity(us,x);
    else
        [P,x] = pdfsample(us,x);
    end

case 'multipdf'
    us=us'; 
    if ischarin('axis',varargin)
        ax = getcharin('axis',varargin);
        varargin = delcharin('axis',varargin);
        x = linspace(ax(1),ax(2),n);
        y = linspace(ax(3),ax(4),n);
%Pout = find(us(:,1)>ax(2) | us(:,1)<ax(1) | us(:,2)>ax(4) | us(:,2)<ax(3));
%Pout = numel(Pout)/size(us,1);
    else
        x = linspace(min(us(:,1))-10*eps,max(us(:,1))+10*eps,n);
        y = linspace(min(us(:,2))-10*eps,max(us(:,2))+10*eps,n);
    end
    [P,x,y] = pdfsample(us,x,y);

otherwise
    error('pas programme')
end

if strcmp(courbe,'multipdf')
    surf(x,y,P+max(max(P))*1.1e-3,varargin{:},'EdgeColor','none','FaceColor','interp')    
    xlim([min(min(x)),max(max(x))]);
    ylim([min(min(y)),max(max(y))]);
%view(-10,87)
elseif tracebar
    bar(x,P,varargin{:})    
    xlim([min(min(x)),max(max(x))]);
    xb = xlim;dx = xb(2)-xb(1);
    xlim([xb(1)-dx/length(x)*2,xb(2)+dx/length(x)*2])
elseif ischarin('fill',varargin)
    plotenveloppe(x,P,varargin{:})   
    xlim([min(min(x)),max(max(x))]);
else
    plot(x,P,varargin{:})
    xlim([min(min(x)),max(max(x))]);
end


if nargout>=1
    varargout{1}=P;
end
if nargout>=2
    varargout{2}=x;
end
if nargout>=3
    varargout{3}=y;
end
