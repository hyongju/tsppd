function [ f,gradf ] = quadcost( x,Q,real_length )
%This function should be used as a function handle for cost function:
%                     f = x'*Q*x
%where in our case Q should be symmetric and x is a vector that store our
%desired permutation matrix.
%   input:  x     --  variable of the function handle.
%           Q     --  The symmetric matrix in f.
%   output: f     --  cost function.
%           gradf --  gradient of f
%   syntax: fun = @(x)quadcost(x,Q,real_length)
    if ~exist('real_length', 'var') || isempty(real_length)
        real_length=length(x);
    end
    f = x(1:real_length)'*Q*x(1:real_length);
    if nargout>1,
       gradf = 2*Q*x(1:real_length);
       gradf = [gradf;zeros(length(x)-real_length,1)];
    end

end
