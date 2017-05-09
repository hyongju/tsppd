function [ f,gradf ] = eig_constr_gene_fun( x,real_length,lambda )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    x_real = x(1:real_length);
    eig_vec = x(real_length+1:end);
    M = [1,x_real'.^2;x_real.^2,x_real*x_real'];
    f = M*eig_vec - lambda*eig_vec;
    gradf = [2*x_real.*eig_vec(2:end);1-lambda;x_real.^2];
    for i = 1:real_length,
       gradf(:,end+1) = [x_real(i)*eig_vec(2:end);x_real(i)^2;x_real(i)*x_real];
       gradf(i,end) = gradf(i,end)+sum(x_real.*eig_vec(2:end))+2*eig_vec(1)*x_real(i);
       gradf(real_length+1+i,end) = gradf(real_length+1+i,end)-lambda;
    end

end

