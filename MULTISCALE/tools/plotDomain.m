function plotDomain(D,varargin)
% function plotDomain(D,D_patch,varargin)
% Display all domains
% D: GEOMOBJECT or MODEL
% D_patch: cell of GEOMOBJECT

p = ImprovedInputParser;
addParameter(p,'legend',true,@islogical);
addParameter(p,'solid',false,@islogical);
addParameter(p,'surface',false,@islogical);
addParameter(p,'view',[],@isnumeric);
addParameter(p,'camup','auto',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
addParameter(p,'Interpreter','latex',@ischar);
if nargin==1 || isempty(varargin) || ischar(varargin{1})
    parse(p,varargin{:})
else
    parse(p,varargin{2:end})
end

varargin = delcharin({'legend','solid','surface','view','camup','FontSize','Interpreter'},varargin);

if ~iscell(D)
    D = {D};
end

if nargin==1 || isempty(varargin) || ischar(varargin{1})
    figure('Name','Domain')
    % set(gcf,'Name','Domain')
    clf
    for l=1:length(D)
        if isa(D{l},'POINT')
            h = plot(D{l},'.');
        elseif isa(D{l},'GEOMOBJECT')
            if p.Results.solid
                h = plot(D{l},'FaceColor',getfacecolor(1),'solid',varargin{:});
            else
                h = plot(D{l},'FaceColor',getfacecolor(1),varargin{:});
            end
        elseif isa(D{l},'MODEL')
            if p.Results.solid
                h = plot(D{l},'FaceColor',getfacecolor(1),'EdgeColor','none','solid',varargin{:});
            else
                h = plot(D{l},'FaceColor',getfacecolor(1),'EdgeColor','none',varargin{:});
            end
        end
        hg = hggroup;
        set(h(:),'Parent',hg);
    end
    if p.Results.legend
        legend(hg,'$\Omega$','Location','NorthEastOutside','Interpreter',p.Results.Interpreter)
    end
    set(gca,'FontSize',p.Results.FontSize)
else
    D_patch = varargin{1};
    if ~iscell(D_patch)
        D_patch = {D_patch};
    end
    n = numel(D_patch);
    figure('Name',['Complement subdomain and patches #' num2str(1:n)])
    % set(gcf,'Name',['Complement subdomain and patches #' num2str(1:n)])
    clf
    hg = cell(1,1+n);
    leg = cell(1,1+n);
    for l=1:length(D)
        if isa(D{l},'POINT')
            h = plot(D{l},'.');
        elseif isa(D{l},'GEOMOBJECT')
            if p.Results.solid
                h = plot(D{l},'FaceColor',getfacecolor(1),'solid',varargin{2:end});
            else
                h = plot(D{l},'FaceColor',getfacecolor(1),varargin{2:end});
            end
        elseif isa(D{l},'MODEL')
            plot(create_boundary(D{l}));
            if p.Results.solid
                h = plot(D{l},'FaceColor',getfacecolor(1),'EdgeColor','none','solid',varargin{2:end});
            else
                h = plot(D{l},'FaceColor',getfacecolor(1),'EdgeColor','none',varargin{2:end});
            end
        end
        hg{1} = hggroup;
        set(h(:),'Parent',hg{1});
    end
    leg{1} = '$\Omega \setminus \Lambda$';
    for k=1:n
        if iscell(D_patch{k})
            for l=1:length(D_patch{k})
                if isa(D_patch{k}{l},'POINT')
                    h = plot(D_patch{k}{l},'.');
                elseif isa(D_patch{k}{l},'GEOMOBJECT')
                    if p.Results.solid
                        h = plot(D_patch{k}{l},'FaceColor',getfacecolor(k+1),'solid',varargin{2:end});
                    else
                        h = plot(D_patch{k}{l},'FaceColor',getfacecolor(k+1),varargin{2:end});
                    end
                elseif isa(D_patch{k}{l},'MODEL')
                    plot(create_boundary(D_patch{k}{l}));
                    if p.Results.solid
                        h = plot(D_patch{k}{l},'FaceColor',getfacecolor(k+1),'EdgeColor','none','solid',varargin{2:end});
                    else
                        h = plot(D_patch{k}{l},'FaceColor',getfacecolor(k+1),'EdgeColor','none',varargin{2:end});
                    end
                end
            end
        else
            if isa(D_patch{k},'POINT')
                h = plot(D_patch{k},'.');
            elseif isa(D_patch{k},'GEOMOBJECT')
                if p.Results.solid
                    h = plot(D_patch{k},'FaceColor',getfacecolor(k+1),'solid',varargin{2:end});
                else
                    h = plot(D_patch{k},'FaceColor',getfacecolor(k+1),varargin{2:end});
                end
            elseif isa(D_patch{k},'MODEL')
                plot(create_boundary(D_patch{k}));
                if p.Results.solid
                    h = plot(D_patch{k},'FaceColor',getfacecolor(k+1),'EdgeColor','none','solid',varargin{2:end});
                else
                    h = plot(D_patch{k},'FaceColor',getfacecolor(k+1),'EdgeColor','none',varargin{2:end});
                end
            end
        end
        hg{k+1} = hggroup;
        set(h(:),'Parent',hg{k+1});
        if n==1
            leg{k+1} = '$\Lambda$';
        else
            leg{k+1} = ['$\Lambda_{' num2str(k) '}$'];
        end
    end
    if p.Results.legend
        legend([hg{:}],leg{:},'Location','NorthEastOutside','Interpreter',p.Results.Interpreter)
    end
    set(gca,'FontSize',p.Results.FontSize)
end

axis image
axis off

if ~isempty(p.Results.view)
    view(p.Results.view)
elseif getindim(D{1})==3 || p.Results.surface
    view(3)
end
if ~isempty(p.Results.camup)
    camup(p.Results.camup)
end
    
end
