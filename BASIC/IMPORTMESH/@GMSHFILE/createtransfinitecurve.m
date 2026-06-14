function u = createtransfinitecurve(u,number,value,varargin)
% function u = createtransfinitecurve(u,number,value,[,meshtype,coeff])

if nargin<2 || isempty(number)
    number = ':';
end

% Flags for options
progression = getcharin('Progression',varargin,[]);
bump = getcharin('Bump',varargin,[]);

options = '';
if ~isempty(progression)
    options = [options ' Using Progression ' tag2str(progression)];
elseif ~isempty(bump)
    options = [options ' Using Bump ' tag2str(bump)];
end

u = stepcounter(u);
s = ['Transfinite Curve ' tagsintobraces(number) ' = ' tag2str(value) options ';\n'];
u = addstring(u,s);
