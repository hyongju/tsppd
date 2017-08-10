function [ y,grady ] = object( x,v,Q,k )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin<4,
    k = 100000;
end
y = x'*Q*x + k*ones(1,v^2)*((x.*(x-1)).^2);
if nargout > 1
    grady = 2*Q*x + k*(4*x.^3-6*x.^2+2*x);
end

end

