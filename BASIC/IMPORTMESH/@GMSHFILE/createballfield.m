function u = createballfield(u,VIn,VOut,XCenter,YCenter,ZCenter,Radius,Thickness,number)
% function u = createballfield(u,VIn,VOut,XCenter,YCenter,ZCenter,Radius,Thickness,number)

if nargin<2 || isempty(VIn)
    VIn = 1e+22;
end
if nargin<3 || isempty(VOut)
    VOut = 1e+22;
end
if nargin<4 || isempty(XCenter)
    XCenter = 0;
end
if nargin<5 || isempty(YCenter)
    YCenter = 0;
end
if nargin<6 || isempty(ZCenter)
    ZCenter = 0;
end
if nargin<7 || isempty(Radius)
    Radius = 0;
end
if nargin<8 || isempty(Thickness)
    Thickness = 0;
end
if nargin<9 || isempty(number)
    number = 1;
end

u = createfield(u,'Box',number);
u = setfieldattribute(u,'VIn',VIn,number);
u = setfieldattribute(u,'VOut',VOut,number);
u = setfieldattribute(u,'XCenter',XCenter,number);
u = setfieldattribute(u,'YCenter',YCenter,number);
u = setfieldattribute(u,'ZCenter',ZCenter,number);
u = setfieldattribute(u,'Radius',Radius,number);
u = setfieldattribute(u,'Thickness',Thickness,number);
