function u = delete(u,name,number,varargin)
% function u = delete(u,name,number[,'recursive'])
% function u = delete(u,name,number[,'recursive'])

% Recursive delete: delete the entities as well as all its sub-entities of lower dimension
recursive = ischarin('recursive',varargin);
% varargin = delonlycharin('recursive',varargin);

if nargin<3 || isempty(number)
    number = ':';
end

% Output capture
prefix = '';
if recursive
    prefix = 'Recursive ';
end

u = stepcounter(u);

entities = [name tagsintobraces(number) ';'];

s = [prefix 'Delete{ ' entities ' }\n'];

u = addstring(u,s);
