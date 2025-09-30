function s = entity(name,number,values)
% function s = entity(name,number,values)

s = [name '(' tag2str(number) ') = ' tagsintobraces(values) ';\n'];
