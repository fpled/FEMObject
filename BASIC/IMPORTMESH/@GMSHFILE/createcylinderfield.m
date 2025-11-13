function u = createcylinderfield(u,VIn,VOut,XCenter,YCenter,ZCenter,XAxis,YAxis,ZAxis,Radius,number)
% function u = createcylinderfield(u,VIn,VOut,XCenter,YCenter,ZCenter,XAxis,YAxis,ZAxis,Radius,number)

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
if nargin<7 || isempty(XAxis)
    XAxis = 0;
end
if nargin<8 || isempty(YAxis)
    YAxis = 0;
end
if nargin<9 || isempty(ZAxis)
    ZAxis = 0;
end
if nargin<10 || isempty(Radius)
    Radius = 0;
end
if nargin<11 || isempty(number)
    number = 1;
end

u = createfield(u,'Box',number);
u = setfieldattribute(u,'VIn',VIn,number);
u = setfieldattribute(u,'VOut',VOut,number);
u = setfieldattribute(u,'XCenter',XCenter,number);
u = setfieldattribute(u,'YCenter',YCenter,number);
u = setfieldattribute(u,'ZCenter',ZCenter,number);
u = setfieldattribute(u,'XAxis',XAxis,number);
u = setfieldattribute(u,'YAxis',YAxis,number);
u = setfieldattribute(u,'ZAxis',ZAxis,number);
u = setfieldattribute(u,'Radius',Radius,number);
