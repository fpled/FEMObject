function [B,detJ] = calc_B(elem,xnode,xgauss)
% function [B,detJ] = calc_B(elem,xnode,xgauss)

nbelem = getnbelem(elem);
nbddl = 3*6;

[detJ,J,Ji,xnode,DNllocal] = calc_detJ(elem,xnode,xgauss);
DNqlocal = getDNq(elem,xgauss);
Nl = getN(elem,xgauss);
DNl = Ji*DNllocal;
DNq = Ji*DNqlocal;

Dxy = [xnode(2,:)-xnode(1,:);xnode(3,:)-xnode(2,:);xnode(1,:)-xnode(3,:)];

% MEMBRANE (UX,UY,RZ)
Zmu = zerosND([3,3*3,sizeND(detJ)]);
Zmv = zerosND([3,3*3,sizeND(detJ)]);
for i=1:3
    Zmu(i,[3*(i-1)+1,mod(3*i,9)+1]) = 1/2;
    Zmu(i,3*(i-1)+3) = -Dxy(i,2)./4;
    Zmu(i,mod(3*i,9)+3) = Dxy(i,2)./4;
    Zmv(i,[3*(i-1)+2,mod(3*i,9)+2]) = 1/2;
    Zmv(i,3*(i-1)+3) = Dxy(i,1)./4;
    Zmv(i,mod(3*i,9)+3) = -Dxy(i,1)./4;
end

DNZmu = DNq(:,4:6)*Zmu;
DNZmv = DNq(:,4:6)*Zmv;

Bm = zerosND([3,3*3,sizeND(detJ)]);
Bm(1,1:3:end) = DNq(1,1:3);
Bm(2,2:3:end) = DNq(2,1:3);
Bm(3,1:3:end) = DNq(2,1:3);
Bm(3,2:3:end) = DNq(1,1:3);
Bm(1,:) = Bm(1,:)+DNZmu(1,:);
Bm(2,:) = Bm(2,:)+DNZmv(2,:);
Bm(3,:) = Bm(3,:)+DNZmu(2,:)+DNZmv(1,:);

% Zmu = zerosND([3,3*3,sizeND(detJ)]);
% Zmv = zerosND([3,3*3,sizeND(detJ)]);
% for i=1:3
%     Zmu(i,3*(i-1)+3) = -Dxy(i,2)./4;
%     Zmu(i,mod(3*i,9)+3) = Dxy(i,2)./4;
%     Zmv(i,3*(i-1)+3) = Dxy(i,1)./4;
%     Zmv(i,mod(3*i,9)+3) = -Dxy(i,1)./4;
% end
% 
% DNZmu = DNq(:,4:6)*Zmu;
% DNZmv = DNq(:,4:6)*Zmv;
% 
% Bm = zerosND([3,3*3,sizeND(detJ)]);
% Bm(1,1:3:end) = DNl(1,:);
% Bm(2,2:3:end) = DNl(2,:);
% Bm(3,1:3:end) = DNl(2,:);
% Bm(3,2:3:end) = DNl(1,:);
% Bm(1,:) = Bm(1,:)+DNZmu(1,:);
% Bm(2,:) = Bm(2,:)+DNZmv(2,:);
% Bm(3,:) = Bm(3,:)+DNZmu(2,:)+DNZmv(1,:);

% DRILLING (UX,UY,RZ)
Bd = zerosND([1,3*3,sizeND(detJ)]);
Bd(1,2:3:end) = DNq(1,1:3)./2;
Bd(1,1:3:end) = -DNq(2,1:3)./2;
Bd = Bd + DNZmv(1,:)./2;
Bd = Bd - DNZmu(2,:)./2;
Bd(1,3:3:end) = - Nl;

% BENDING (UZ,RX,RY)=(UZ,-BY,BX)
l2 = (Dxy(:,1).^2 + Dxy(:,2).^2);
l = sqrt(l2);

a = -3/2*Dxy(:,1)./l2;
b = (1/4*Dxy(:,1).^2-1/2*Dxy(:,2).^2)./l2;
c = 3/4*(Dxy(:,1).*Dxy(:,2))./l2;
d = -3/2*Dxy(:,2)./l2;
e = (1/4*Dxy(:,2).^2-1/2*Dxy(:,1).^2)./l2;

Zbx = zerosND([3,3*3,sizeND(detJ)]);
Zby = zerosND([3,3*3,sizeND(detJ)]);
for i=1:3
    Zbx(i,[3*(i-1)+[1:3],mod(3*i,9)+[1:3]]) = -[a(i),b(i),c(i),-a(i),b(i),c(i)];
    Zby(i,[3*(i-1)+[1:3],mod(3*i,9)+[1:3]]) = -[d(i),c(i),e(i),-d(i),c(i),e(i)];
end

DNZbx = DNq(:,4:6)*Zbx;
DNZby = DNq(:,4:6)*Zby;

Bb = zerosND([3,3*3,sizeND(detJ)]);
Bb(1,2:3:end) = DNq(1,1:3);
Bb(2,3:3:end) = DNq(2,1:3);
Bb(3,2:3:end) = DNq(2,1:3);
Bb(3,3:3:end) = DNq(1,1:3);
Bb(1,:) = Bb(1,:)+DNZbx(1,:);
Bb(2,:) = Bb(2,:)+DNZby(2,:);
Bb(3,:) = Bb(3,:)+DNZbx(2,:)+DNZby(1,:);

% a = -3/2*Dxy(:,1)./l2;
% b = 3/4*Dxy(:,1).^2./l2;
% c = 3/4*(Dxy(:,1).*Dxy(:,2))./l2;
% d = -3/2*Dxy(:,2)./l2;
% e = 3/4*Dxy(:,2).^2./l2;
% 
% Zbx = zerosND([3,3*3,sizeND(detJ)]);
% Zby = zerosND([3,3*3,sizeND(detJ)]);
% for i=1:3
%     Zbx(i,[3*(i-1)+[1:3],mod(3*i,9)+[1:3]]) = -[a(i),b(i),c(i),-a(i),b(i),c(i)];
%     Zby(i,[3*(i-1)+[1:3],mod(3*i,9)+[1:3]]) = -[d(i),c(i),e(i),-d(i),c(i),e(i)];
% end
% 
% DNZbx = DNq(:,4:6)*Zbx;
% DNZby = DNq(:,4:6)*Zby;
% 
% Bb = zerosND([3,3*3,sizeND(detJ)]);
% Bb(1,2:3:end) = DNl(1,:);
% Bb(2,3:3:end) = DNl(2,:);
% Bb(3,2:3:end) = DNl(2,:);
% Bb(3,3:3:end) = DNl(1,:);
% Bb(1,:) = Bb(1,:)+DNZbx(1,:);
% Bb(2,:) = Bb(2,:)+DNZby(2,:);
% Bb(3,:) = Bb(3,:)+DNZbx(2,:)+DNZby(1,:);

Bbtemp = Bb;
Bb(:,2:3:end) = -Bbtemp(:,3:3:end); % RX = -BY
Bb(:,3:3:end) = Bbtemp(:,2:3:end); % RY = BX

% TOTAL (UX,UY,UZ,RX,RY,RZ)
B = zerosND([7,3*6,sizeND(detJ)]);
repm=[];repb=[];
for i=1:3
    repm = [repm,6*(i-1)+[1,2,6]];
    repb = [repb,6*(i-1)+[3,4,5]];
end
B(1:3,repm) = Bm;
B(4:6,repb) = Bb;
B(7,repm) = Bd;

% ROTATION
P = calc_P(elem);
B = B*P;

if nargout>2
    x = calc_x(elem,xnode,xgauss);
end
