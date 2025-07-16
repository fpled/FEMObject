function varargout = plot(u,varargin)
% function [P,x] = plot(u,courbe)
% trace la pdf ou cdf ou icdf de u (PCMATRIX)
% courbe  = 'pdf' ou 'icdf' ou 'cdf' ou 'fun' ('fun' par defaut)
% 'fun' trace la surface de reponse en fonction des variables de base du chaos
% plot(u,courbe,'npts',n) trace la courbe avec n points (100 par defaut)
% plot(u,courbe,'bar','npts',n) trace des barres.
% plot(u,courbe,'nbs',N) nombre de tirage (10000 par defaut)
% tous les autres arguments sont passes en argument de plot ou bar
%
% See also sampleplot, PCMATRIX/pdfplot, PCMATRIX/cdfplot, PCMATRIX/funplot, RANDVAR/plot

if any(ischarin({'pdf','cdf','icdf','multipdf'},varargin))
    N=getcharin('nbs',varargin,1e4);
    varargin = delcharin({'nbs'},varargin);
    us=randomlimit(u,N);
    
    if numel(u)==2 & ischarin('multipdf',varargin)
        us = us;
    elseif numel(u)~=1
        error('tracer une multipdf pour 2 variables et pas plus')
    else
        us = us(:);
    end
    
    [P,x]=sampleplot(us,varargin{:});
else
    RV = RANDVARS(u);
    n=getcharin('npts',varargin,100);
    
    switch length(RV)
        case 1
            if ischarin('axis',varargin)
                x = getcharin('axis',varargin);
                varargin=delcharin('axis',varargin);
            else
                x = getdomainborne(RV{1});
            end
            x = linspace(x(1),x(2),n);
            P = randomeval(u,x(:));
            varargin=delonlycharin('bar',varargin);
            varargin = delcharin({'npts','nbs','nbbar'},varargin);
            plot(x,P,varargin{:})
            
        case 2
            if ischarin('axis',varargin)
                ax = getcharin('axis',varargin);;
                x = ax(1:2);
                y = ax(3:4);
                varargin=delcharin('axis',varargin);
            else
                x = getdomainborne(RV{1});
                y = getdomainborne(RV{2});
            end
            ax = [x(1),x(2),y(1),y(2)]
            x = linspace(x(1),x(2),n);
            y = linspace(y(1),y(2),n);
            [x,y] = meshgrid(x,y);
            P = randomeval(u,[x(:),y(:)]);
            P = reshape(P,n,n);
            surf(x,y,P,'EdgeColor','none','FaceColor','interp',varargin{:})
            axz = [min(min(P)),max(max(P))];
            axis([ax,axz])
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
end

