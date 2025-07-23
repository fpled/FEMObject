function M = MODEL(mode)
%function M = MODEL(mode)
% constructeur de la classe MODEL
% mode : UNID, PLANE, TRID
if nargin==0
    M = MODEL(2);
elseif nargin==1 && isa(mode,'MODEL')
    M = mode;
else
    if isa(mode,'char')
        M.mode = mode;
        M.dim = getindim(mode);
    elseif isa(mode,'double')
        indim = mode;
        M.mode = getmode(indim);
        M.dim = indim;
    end
    
    switch M.mode
        case 'UNID'
            M.syscoord = CARTESIAN1D();
        case 'PLANE'
            M.syscoord = CARTESIAN2D();
        case 'TRID'
            M.syscoord = CARTESIAN3D();
    end
    
    
    M.nbelem = 0;
    M.nbgroupelem = 0;
    M.groupelem = cell(0,1);
    M.repelemingroupelem = zeros(0,2);
    
    M.nbnode = 0 ;
    M.node = NODE(M.syscoord);
    M.connec = struct();
    
    M.nbddl = 0;
    M.facets = cell(1,0);
    M.ridges = cell(1,0);
    M.peaks = cell(1,0);
    
    
    M=class(M,'MODEL',BCOND(M.node));
    
end


function indim = getindim(mode)

switch mode
    case 'UNID'
        indim = 1;
    case 'PLANE'
        indim = 2;
    case 'TRID'
        indim = 3;
end

function mode = getmode(indim)

switch indim
    case 1
        mode = 'UNID';
    case 2
        mode = 'PLANE';
    case 3
        mode = 'TRID';
end

