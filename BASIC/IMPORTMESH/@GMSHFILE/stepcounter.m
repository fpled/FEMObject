function u = stepcounter(u,number)
% function u = stepcounter(u,number)

if nargin==1 || ischar(number)
    u.counter = u.counter+1;
else
    u.counter = max(number,u.counter+1);
end
