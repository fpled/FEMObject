function w = expand3D(u)
% function w = expand3D(u)

s = num2cell(size(u));
if s{1}==2
    u3 = zerosND(1,s{2:end});
    u.MYDOUBLEND = [u.MYDOUBLEND;u3];
end

w = u;
