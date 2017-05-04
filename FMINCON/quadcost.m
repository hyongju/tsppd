function [ f,gradf ] = quadcost( x,Q )
%This function should be used as a function handle for cost function:
%                     f = x'*Q*x
%where in our case Q should be symmetric and x is a vector that store our
%desired permutation matrix.
%   input:  x     --  variable of the function handle.
%           Q     --  The symmetric matrix in f.
%   output: f     --  cost function.
%           gradf --  gradient of f
%   syntax: fun = @(x)quadcost(x,Q)
    f = x'*Q*x;
    if nargout>1,
       gradf = 2*Q*x; 
    end

end
