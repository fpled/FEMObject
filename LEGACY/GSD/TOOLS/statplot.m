function varargout = plot(u,varargin)
% function [P,x,y] = plot(u,courbe)
% trace la pdf ou cdf ou icdf de u
% courbe  = 'pdf' ou 'icdf' ou 'cdf' ou 'fun' ('fun' par defaut)
% 'fun' trace la surface de reponse en fonction des variables de base du chaos
% plot(u,courbe,'npts',n) trace la courbe avec n points (100 par defaut)
% plot(u,courbe,'bar','npts',n) trace des barres. 
% plot(u,courbe,'nbs',N) nombre de tirage (10000 par defaut)    
% tous les autres arguments sont passes en argument de plot ou bar
%
% See also sampleplot, pdfplot, cdfplot, funplot, plot


if any(ischarin({'pdf','cdf','icdf','multipdf'},varargin)) %% COURBES STATS

    if israndom(u)
        N=getcharin('nbs',varargin,1e4);
        varargin = delcharin({'nbs'},varargin);
        us=randomlimit(u,N);
        num = numel(u);
    else
        us=u;
    end

    if num==2 && ischarin('multipdf',varargin)  
%ok
    elseif numel(u)~=1 
        error('tracer une multipdf pour 2 variables et pas plus')
    else
        us = us(:);
    end
    varargout=cell(1,nargout);
    [varargout{:}]=sampleplot(us,varargin{:});


else   %  SURFACE DE REPONSE

    if ~israndom(u)
        error('surface de reponse pas possible avec des echantillons')
    end

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

        if isclassin('RANDVARS',varargin);
            RVplot = getclassin('RANDVARS',varargin);
            varargin = delclassin('RANDVARS',varargin);
            x = transfer(RV,RVplot,x(:));
            plot(x,P,varargin{:})
        else
            plot(x,P,varargin{:})
        end
    case 2
        varargin = delcharin({'npts','nbs','nbbar'},varargin);        
        if ischarin('axis',varargin)
            ax = getcharin('axis',varargin);
            x = ax(1:2);
            y = ax(3:4);
            varargin=delcharin('axis',varargin);
        else
            x = getdomainborne(RV{1});
            y = getdomainborne(RV{2});
        end
        ax = [x(1),x(2),y(1),y(2)];
        x = linspace(x(1),x(2),n);
        y = linspace(y(1),y(2),n);
        [x,y] = meshgrid(x,y);
        P = randomeval(u,[x(:),y(:)]);
        P = reshape(P,n,n);
        if isclassin('RANDVARS',varargin);
            RVplot = getclassin('RANDVARS',varargin);
            xplot = getdomainborne(RVplot{1});
            yplot = getdomainborne(RVplot{2});
            ax = [xplot(1),xplot(2),yplot(1),yplot(2)];
            varargin = delclassin('RANDVARS',varargin);
            xy = transfer(RV,RVplot,[x(:),y(:)]);
            x=reshape(xy(:,1),size(x));
            y=reshape(xy(:,2),size(y));
            if ~ischarin('FaceColor',varargin)
                varargin = setcharin('FaceColor',varargin,'interp');
            end
            if ~ischarin('EdgeColor',varargin)
                varargin = setcharin('EdgeColor',varargin,'none');
            end
            surf(x,y,P,varargin{:})    
        else
            if ~ischarin('FaceColor',varargin)
                varargin = setcharin('FaceColor',varargin,'interp');
            end
            if ~ischarin('EdgeColor',varargin)
                varargin = setcharin('EdgeColor',varargin,'none');
            end
            surf(x,y,P,varargin{:})    
        end
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

