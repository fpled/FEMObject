function mat = FOUR_ISOT(varargin)
% function mat = FOUR_ISOT('k',k,'c',c,'k2',k2,'b',b,'r',r,'r2',r2,'r3',r3)
% equation d'advection-diffusion-reaction transitoire non-lineaire
% c du/dt + div(q) + b grad(u) + r*u + r2*u^2 + r3*u^3 = f
% relation entre le flux q et la temperature u
% q = -(k+k2*u^2) grad(u)
% 'k' coefficient de diffusion (scalaire) intervenant dans le terme de diffusion lineaire : -k grad(u)
% 'k2' coefficient de diffusion (scalaire) intervenant dans le terme de diffusion non-lineaire : -k2*u^2 grad(u)
% 'b' vitesse de transport (vecteur) intervenant dans le terme d'advection : b grad(u)
% 'c' conductivite (scalaire) intervenant dans le terme d'evolution : c du/dt
% 'r' parametre de reaction (scalaire) intervenant dans le terme de reaction lineaire : r*u
% 'r2' parametre de reaction (scalaire) intervenant dans le terme de reaction non-lineaire : r*u^2
% 'r3' parametre de reaction (scalaire) intervenant dans le terme de reaction non-lineaire : r*u^3

mat = struct();

param = struct(varargin{:});

mat.k2 = ischarin('k2',varargin);
mat.b  = ischarin('b',varargin);
mat.r  = ischarin('r',varargin);
mat.r2 = ischarin('r2',varargin);
mat.r3 = ischarin('r3',varargin);

if ~ischarin('stabilize',varargin)
    param.stabilize = 0;
end

if ~ischarin('PFregularization',varargin)
    param.PFregularization = 'AT2';
end

if isfield(param,'lcorr')
    if ~isfield(param,'shinozuka')
        param.shinozuka = @(x) shinozukaFun(x);
    end
end

matp = MATERIAL('FOUR_ISOT',param);
mat = class(mat,'FOUR_ISOT',matp);

