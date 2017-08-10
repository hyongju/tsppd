function [ A ] = mat_col_switch( A,idx1,idx2 )
%this function helps to switch two columns in a matrix A
%   Detailed explanation goes here
if max(idx1,idx2)>size(A,2),
   error('the maximum index exceeds the number of columns in A');
end
temp = A(:,idx1);
A(:,idx1) = A(:,idx2);
A(:,idx2) = temp;


end

