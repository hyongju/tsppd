function [ f_hess ] = ehess(u,v,Q,A,Base1,Base2,P,N,aux)
%This function generates the Euclidian Hessian of the cost function. (user
%design)
%   Input: user indecates
    n = length(Base1);
    f_hess = 0;
    for i = 1:n,
        for j = 1:n,
            AA = Base2{j}'*A'*A*Base2{i};
            BB = Base1{i}*Base1{j}';
            f_hess = f_hess + AA*v*u'*BB*u + BB*u*v'*AA*u + AA*u*u'*BB*v;
            f_hess = f_hess + BB'*v*u'*AA'*u + AA'*u*v'*BB'*u + BB'*u*u'*AA'*v;
        end
    end
    f_hess = 2*N*f_hess + 2*aux*Q*aux'*v - 2*N*P*v; 
    

end

