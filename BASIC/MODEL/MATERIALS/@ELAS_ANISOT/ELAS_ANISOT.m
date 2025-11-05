function mat = ELAS_ANISOT(varargin)
% function mat = ELAS_ANISOT('parameter',value)
% Parameters:
% 'C' elasticity matrix
% 'RHO' mass density
% 'DIM3' thickness in the 2D case
% 'S' cross-section area in the 1D case
% 'MU' damping factor

mat = struct();

param = struct(varargin{:});

if ~isfield(param,'DIM3')
    param.DIM3 = 1;
end

if isfield(param,'d')
    if ~isfield(param,'PFM')
        param.PFM = 'HeFreddi';
    end
    if ~isfield(param,'PFS')
        param.PFS = 'Strain';
    end
    if ~isfield(param,'g')
        if isfield(param,'h')
            % param.g = @(d,h) (1-d.*(1-h)).^2;
            param.g = @(d,h) (1-d+h).^2;
        else
            param.g = @(d) (1-d).^2;
        end
    end
    if ~isfield(param,'k')
        param.k = 0;
    end
end

if isfield(param,'lcorr')
    if ~isfield(param,'shinozuka')
        param.shinozuka = @(x) shinozukaFun(x);
    end
end

matp = MATERIAL('ELAS_ANISOT',param);
mat = class(mat,'ELAS_ANISOT',matp);

