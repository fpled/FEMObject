function ls=LSHYPERPLANE(dim,varargin)
% function ls=LSHYPERPLANE(dim,Px,Vx)   un point et un vecteur
% function ls=LSHYPERPLANE(dim,Px,Py,Vx,Vy)
% function ls=LSHYPERPLANE(dim,Px,Py,Pz,Vx,Vy,Vz)

ls=struct();
param = cell(6,2);
param(:,1)={'Px','Py','Pz','Vx','Vy','Vz'};

rep = [1:dim,3+[1:dim]]; 

param = param(rep,:);
param(:,2)=varargin(1:length(rep));

lsp = LEVELSET(param,varargin{:});
ls=class(ls,'LSHYPERPLANE',lsp);
