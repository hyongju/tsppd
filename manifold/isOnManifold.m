function [ res ] = isOnManifold( y,A,b )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    res = 1;
    x = y(2:end);
    for i = 1:length(A)
       tempA = A{i}(2:end,2:end);
       if abs(tempA*x-1)>=1e-7,
           res = 0;
           break;
       end
    end



end

