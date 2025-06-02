function [y,s] = calc_response(fun,x,n,varargin)
% function [y,s] = calc_response(fun,x,n,varargin)
% Compute response matrix y returned by function fun applied to sample matrix x
% fun : function_handle of random variables
% x : sample matrix
% n : number of output arguments for fun
% varargin : arguments of function fun
% call y = fun(x,varargin{:}) for sample matrix x
% x(:,i) is the sample vector containing all the evaluations of random variable i
% x(k,:) is the sample vector containing the k-th evaluation of all random variables
% y is the response matrix containing the values returned by function fun applied to sample matrix x

if ischarin('display',varargin)
    varargin = delonlycharin('display',varargin);
end

if nargin<3 || isempty(n)
    n = 1;
end

y = cell(1,n);
s = cell(1,n);
[y{:}] = fun(x,varargin{:});
for i=1:n
    if isa(y{i},'cell')
        for l=1:length(y{i})
            if isa(y{i}{l},'MULTIMATRIX')
                s{i}{l} = size(y{i}{l});
                y{i}{l} = double(y{i}{l});
            else
                s{i}{l} = [size(y{i}{l},1),1];
            end
        end
    else
        if isa(y{i},'MULTIMATRIX')
            s{i} = size(y{i});
            y{i} = double(y{i});
        else
            s{i} = [size(y{i},1),1];
        end
    end
end

end
