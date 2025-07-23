function [Pbase,Ptop] = getcenters(C)
% function [Pbase,Ptop] = getcenters(C)

Pbase = getc(C);
Ptop  = getctop(C);

Pbase = POINT(Pbase);
Ptop  = POINT(Ptop);
