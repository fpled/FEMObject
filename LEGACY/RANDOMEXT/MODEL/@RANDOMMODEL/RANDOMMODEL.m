function M = RANDOMMODEL(model)
%function M = RANDOMMODEL(model)
% constructeur de la classe RANDOMMODEL
% mode : UNID, PLANE, TRID
if nargin==0 
    M=struct();
    M=class(M,'RANDOMMODEL',MODEL());
elseif nargin==1 
    if isa(model,'RANDOMMODEL')
        M = model;
    elseif isa(model,'MODEL')
        M=struct();
        M=class(M,'RANDOMMODEL',model);
    end
end

