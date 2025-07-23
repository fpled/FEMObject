function display(D)
% function display(D)

disp(' ')
disp([inputname(1) ' = (' class(D) ')'])
disp(' ')
display(D.P{1})
for k=1:length(D.V)
    display(D.V{k})
end
