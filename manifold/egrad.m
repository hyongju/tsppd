function [ f_grad ] = egrad(u,Q,A,Base1,Base2,P,N,aux)
%This function generates the Euclidian gradient of the cost function. (user
%design)
%   Input: user indecates
    
    f_grad = 0;
    n = length(Base1);
    for i = 1:n,
        for j = 1:n,
            f_grad = f_grad + Base1{j}*Base1{i}'*u*u'*Base2{i}'*A'*A*Base2{j};
        end
    end
    f_grad = 2*N*(f_grad+f_grad')*u + 2*aux*Q*aux'*u-2*N*P*u;

end

