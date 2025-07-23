function S = plus(S,u)
% function S1 = plus(S,u)

if isa(u,'VECTOR')

   S.node = S.node+u;

elseif isa(u,'MYDOUBLE') || isa(u,'double') 
   
   u = unfreevector(S.BCOND,double(u));

   d = findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS')));
   u = u(d);
   u = reshape(u,getindim(S.syscoord),size(u,1)/getindim(S.syscoord));

   S.node = S.node + VECTOR(u);

end