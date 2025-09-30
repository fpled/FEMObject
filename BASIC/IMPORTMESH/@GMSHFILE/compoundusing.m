function u = compoundusing(u,name,number,values,numIntervals)
% function u = compoundusing(u,name,number,values,numIntervals)

u = stepcounter(u);
s = ['Compound ' name '(' tag2str(number) ') = ' tagsintobraces(values) ' Using ' tag2str(numIntervals) ';\n'];
u = addstring(u,s);
