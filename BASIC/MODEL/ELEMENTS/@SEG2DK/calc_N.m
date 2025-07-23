function [N,detJ,x] = calc_N(elem,xnode,xgauss)
% function [N,detJ,x] = calc_N(elem,xnode,xgauss)

nbelem = getnbelem(elem);
nbddl = 2*6;

[detJ,J,Ji] = calc_detJ(elem,xnode,xgauss);

sys = getsyscoord(elem);
base = getbase(sys);
eX = VECTOR(base(:,1));
eY = VECTOR(base(:,2));
sys = CARTESIAN2D(eX,eY);
xnode = xnode*sys;

Nl = MYDOUBLEND(getN(elem,xgauss));
Nq = MYDOUBLEND(getNq(elem,xgauss));

Dxy = xnode(2,:)-xnode(1,:);

% MEMBRANE (UX,UY,RZ)
Zmu = zerosND([1,2*3,sizeND(detJ)]);
Zmv = zerosND([1,2*3,sizeND(detJ)]);
Zmu(1,[1,4]) = 1/2;
Zmu(1,3) = -Dxy(1,2)./4;
Zmu(1,6) = Dxy(1,2)./4;
Zmv(1,[2,5]) = 1/2;
Zmv(1,3) = Dxy(1,1)./4;
Zmv(1,6) = -Dxy(1,1)./4;

NZmu = Nq(:,3)*Zmu;
NZmv = Nq(:,3)*Zmv;

Nm = zerosND([3,2*3,sizeND(detJ)]);
Nm(1,1:3:end) = Nq(1,1:2);
Nm(2,2:3:end) = Nq(1,1:2);
Nm(1,:) = Nm(1,:)+NZmu(1,:);
Nm(2,:) = Nm(2,:)+NZmv(1,:);
Nm(3,3:3:end) = Nl;

% Zmu = zerosND([1,2*3,sizeND(detJ)]);
% Zmv = zerosND([1,2*3,sizeND(detJ)]);
% Zmu(1,3) = -Dxy(1,2)./4;
% Zmu(1,6) = Dxy(1,2)./4;
% Zmv(1,3) = Dxy(1,1)./4;
% Zmv(1,6) = -Dxy(1,1)./4;
% 
% NZmu = Nq(:,3)*Zmu;
% NZmv = Nq(:,3)*Zmv;
% 
% Nm = zerosND([3,2*3,sizeND(detJ)]);
% Nm(1,1:3:end) = Nl;
% Nm(2,2:3:end) = Nl;
% Nm(1,:) = Nm(1,:)+NZmu(1,:);
% Nm(2,:) = Nm(2,:)+NZmv(1,:);
% Nm(3,3:3:end) = Nl;

% BENDING (UZ,RX,RY)=(UZ,-BY,BX)
l2 = (Dxy(:,1).^2 + Dxy(:,2).^2);
l = sqrt(l2);

a = -3/2*Dxy(:,1)./l2;
b = (1/4*Dxy(:,1).^2-1/2*Dxy(:,2).^2)./l2;
c = 3/4*(Dxy(:,1).*Dxy(:,2))./l2;
d = -3/2*Dxy(:,2)./l2;
e = (1/4*Dxy(:,2).^2-1/2*Dxy(:,1).^2)./l2;

Zbx = zerosND([1,2*3,sizeND(detJ)]);
Zby = zerosND([1,2*3,sizeND(detJ)]);
Zbx(1,:) = -[a,b,c,-a,b,c];
Zby(1,:) = -[d,c,e,-d,c,e];

NZbx = Nq(:,3)*Zbx;
NZby = Nq(:,3)*Zby;

Nb = zerosND([3,2*3,sizeND(detJ)]);
Nb(1,1:3:end) = Nl;
Nb(2,2:3:end) = Nq(1,1:2);
Nb(3,3:3:end) = Nq(1,1:2);
Nb(2,:) = Nb(2,:)+NZbx(1,:);
Nb(3,:) = Nb(3,:)+NZby(1,:);

% a = -3/2*Dxy(:,1)./l2;
% b = 3/4*Dxy(:,1).^2./l2;
% c = 3/4*(Dxy(:,1).*Dxy(:,2))./l2;
% d = -3/2*Dxy(:,2)./l2;
% e = 3/4*Dxy(:,2).^2./l2;
% 
% Zbx = zerosND([1,2*3,sizeND(detJ)]);
% Zby = zerosND([1,2*3,sizeND(detJ)]);
% Zbx(1,:) = -[a,b,c,-a,b,c];
% Zby(1,:) = -[d,c,e,-d,c,e];
% 
% NZbx = Nq(:,3)*Zbx;
% NZby = Nq(:,3)*Zby;
% 
% Nb = zerosND([3,2*3,sizeND(detJ)]);
% Nb(1,1:3:end) = Nl;
% Nb(2,2:3:end) = Nl;
% Nb(3,3:3:end) = Nl;
% Nb(2,:) = Nb(2,:)+NZbx(1,:);
% Nb(3,:) = Nb(3,:)+NZby(1,:);

Nbtemp = Nb ;
Nb(:,2:3:end) = -Nbtemp(:,3:3:end); % RX = -BY
Nb(:,3:3:end) = Nbtemp(:,2:3:end); % RY = BX

Nbtemp = Nb;
Nb(2,:) = -Nbtemp(3,:); % RX = -BY
Nb(3,:) = Nbtemp(2,:); % RY = BX

% TOTAL (UX,UY,UZ,RX,RY,RZ)
N = zerosND([6,2*6,sizeND(detJ)]);
repm=[];repb=[];
for k=1:2
    repm = [repm,6*(k-1)+[1,2,6]];
    repb = [repb,6*(k-1)+[3,4,5]];
end
N([1,2,6],repm) = Nm;
N([3,4,5],repb) = Nb;

% ROTATION
P = calc_P(elem);
N = N*P;

if nargout>2
    x = calc_x(elem,xnode,xgauss);
end
