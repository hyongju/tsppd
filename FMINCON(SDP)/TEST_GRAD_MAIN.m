function [ ] = TEST_GRAD_MAIN( cost,constr,v )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    disp('test gradients: START.');
    
    res = grad_test(cost,v^2);
    if res==1,
        disp('test cost funcion: valid.');
    else
        disp('test cost funcion: invalid.');
    end
    
    disp('testing constraints...');
    n = length(constr.eq.A);
    for i = 1:n,
%         i
        for j = 1:size(constr.eq.A{i},1),
%             j
            temp = @(x)generate_fun(x,constr.eq.A{i}(j,:),constr.eq.B{1}(j,:));
            res = grad_test(temp,v^2);
            if res==0,
               disp(['gradient fails for the ',num2str(j),'th equality constraint']); 
            end
        end
    end
    n = length(constr.ineq.A);
    for i = 1:n,
%         i
        for j = 1:size(constr.ineq.A{i},1),
%             j
            temp = @(x)generate_fun(x,constr.ineq.A{i}(j,:),constr.ineq.B{1}(j,:));
            res = grad_test(temp,v^2);
            if res==0,
               disp(['gradient fails for the ',num2str(j),'th inequality constraint']); 
            end
        end
    end
    
    disp('testing eigenvalue constraint');
    temp = @(x)eig_constr_gene_fun(x,v,v+1);
    res = grad_test(temp,v+v+1);
    if res==0,
        disp('gradient fails for eigenvalue constraint');
    end
    disp('finish testing constraints.');
    disp('test gradients: END.');

end

