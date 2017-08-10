function [ f ] = cost( u,Q,A,b,N,aux )
%This function generates a cost function for a manopt problem (user design)
%   Input: user indecates

    X = u*u';
    f = trace(aux'*X*aux*Q) + N*norm(A*reshape(X,[],1)-b,'fro')^2;

end

