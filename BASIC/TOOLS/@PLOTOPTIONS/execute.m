function execute(P,varargin)


if ischarin('after',varargin)
    makepause = getparam(P.PARAMETERS,'pause');
    plotiter = getparam(P.PARAMETERS,'plotiter');
    plottime = getparam(P.PARAMETERS,'plottime');
    xli = getparam(P.PARAMETERS,'xlim');
    yli = getparam(P.PARAMETERS,'ylim');
    zli = getparam(P.PARAMETERS,'zlim');
    noxtick = getparam(P.PARAMETERS,'noxtick');
    noytick = getparam(P.PARAMETERS,'noytick');
    noztick = getparam(P.PARAMETERS,'noztick');
    boxon = getparam(P.PARAMETERS,'boxon');
    boxstylefull = getparam(P.PARAMETERS,'boxstylefull');
    
    pausetime = getparam(P.PARAMETERS,'pausetime');
    colorma = getparam(P.PARAMETERS,'colormap');
    colorba = getparam(P.PARAMETERS,'colorbar');
    if ~isa(colorma,'char') || ~strcmp('colorma','default')
        colmap = true;
    else
        colmap = false;
    end
    ax = getparam(P.PARAMETERS,'axis');
    cax = getparam(P.PARAMETERS,'caxis');
    axison = getparam(P.PARAMETERS,'axison');
    axissquare = getparam(P.PARAMETERS,'axissquare');
    axisimage = getparam(P.PARAMETERS,'axisimage');
    numview = getparam(P.PARAMETERS,'view');
    up_vector = getparam(P.PARAMETERS,'camup');
    camera_position = getparam(P.PARAMETERS,'campos');
    fontsize = getparam(P.PARAMETERS,'FontSize');
    
    if ~isempty(ax)
        axis(ax);
    end
    
    if ~isempty(cax)
        caxis(cax);
    end
    
    if ~isempty(xli)
        xlim(xli)
    end
    if ~isempty(yli)
        ylim(yli)
    end
    if ~isempty(zli)
        zlim(zli)
    end
    if noxtick
        set(gca,'XTick',[])
    end
    if noytick
        set(gca,'YTick',[])
    end
    if noztick
        set(gca,'ZTick',[])
    end
    
    if axison
        axis on
    else
        axis off
    end
    
    if axissquare
        axis square
    elseif axisimage
        axis image
    end
    
    if ~isempty(getparam(P.PARAMETERS,'pbaspect'))
        daspect(getparam(P.PARAMETERS,'pbaspect'));
        pbaspect(getparam(P.PARAMETERS,'pbaspect'));
    end
    
    if ~isempty(getparam(P.PARAMETERS,'daspect'))
        daspect(getparam(P.PARAMETERS,'daspect'));
        pbaspect(getparam(P.PARAMETERS,'daspect'));
    end
    
    if boxon
        box on
    end
    if boxstylefull
        set(gca,'BoxStyle','full')
    end
    
    if ~isempty(numview)
        view(numview)
    end
    if ~isempty(up_vector)
        camup(up_vector)
    end
    if ~isempty(camera_position)
        campos(camera_position);
    end
    if colmap
        colormap(colorma)
    end
    if colorba
        colorbar
    end
    if fontsize
        set(gca,'FontSize',fontsize)
    end
    
    if makepause
        pause(pausetime)
    end
    
end

if ischarin('before',varargin)
    
    clearfig = getparam(P.PARAMETERS,'clf');
    if clearfig
        clf
    end
    
end
