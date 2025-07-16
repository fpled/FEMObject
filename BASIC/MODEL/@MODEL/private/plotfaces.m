function varargout = plotfaces(faces,varargin)

islegend = getcharin('legend',varargin,false);

if nargin==1 || ~isa(varargin{1},'double') || (length(varargin{1})==1 && varargin{1}==0)
    scanfaces = 1:length(faces);
else
    scanfaces=varargin{1};
end

options = varargin;
Handles=[];
leg = cell(1,0);
for i=1:length(scanfaces)
    if getgroupelemdim(faces{scanfaces(i)})==0
        options = setcharin('MarkerEdgeColor',options,getfacecolor(scanfaces(i)));
        options = addcharin('numelem',options);
    else
        if ~ischarin('FaceColor',varargin)
            options = setcharin('FaceColor',options,getfacecolor(scanfaces(i)));
        end
    end
    
    Handlestemp = plot(faces{scanfaces(i)},options{:});
    Handles = [Handles , Handlestemp(1)];
    for j=1:length(Handlestemp(1))
        leg = [leg , {num2str(scanfaces(i))}];
    end
end

if ~isempty(Handles) && islegend
    legend(Handles,leg{:});
end

if nargout>=1
    varargout{1}=Handles;
end
