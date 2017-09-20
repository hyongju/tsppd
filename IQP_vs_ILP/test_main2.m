%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET INITIAL GUESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc


n = 6 ;             % number of custumers(n)
k = 1;              % capacity   
v = 2*n+1;
beta = 0.5;
appx1 = 0.1;
appx2 = 1e-1;
% alpha = 0.01;
% construct map and cost map
% file = strcat('map',num2str(n),'.mat');   
% if ~exist(file),




nIter = 1;

for m = 1:nIter
    m
    rng('shuffle');     % random seed: shuffle
    vert = rand(2*n+1,2);
    % else
    %     load(file);
    % end

    % cost matrix
    c = zeros(v,v);
    for i = 1:size(c,1)
        for j = 1:size(c,2)
            c(i,j) = norm(vert(i,:)-vert(j,:));
        end
    end
    % c = c + eye(v); % add identity to prevent zeros from appearing in the diagonal

    % make c PSD...

    for i = 1:size(c,1)
        sum_c(i) = 0;
        for j = 1:size(c,2)
            sum_c(i) = sum_c(i) + abs(c(i,j));
        end
        c(i,i) = c(i,i) + sum_c(i);
        cM(i) = c(i,i);
    end
    for i = 1:size(c,1)
        c(i,i) = beta*max(cM);
    end
    all(eig(c)>=0)  % check if c is PSD



    [init_ILP,init_IQP,obj_IG] = calcInitalValue(n,k,v,c,appx1);

    [outputILP,obj_ILP,solution_ILP,outputIQP,obj_IQP,solution_IQP] = main_hp_fu6nc2(n,k,v,c,appx2,init_ILP,init_IQP);


    
    
    
%     [outputILP,obj_ILP,solution_ILP,outputIQP,obj_IQP,solution_IQP] = main_hp_func(n,k,i,appx);
    timeCmp(m,1) = outputILP.solvertime;
    timeCmp(m,2) = outputIQP.solvertime;
    timeCmp(m,3) = obj_ILP;
    timeCmp(m,4) = obj_IQP;
    timeCmp(m,5) = obj_IG;
end


figure,
plot(1:nIter,timeCmp(:,1),'bo-',1:nIter,timeCmp(:,2),'rs-');legend('ILP','IQP');
hold on;
title(sprintf('n = %d, k = %d', n,k));
set(gca,'FontSize',16);
xlabel('sample #');ylabel('solver time');


figure,
hist(timeCmp(:,1:2));hold on;
title(sprintf('n = %d, k = %d', n,k));
legend('ILP','IQP');
xlabel('solver time');
ylabel('number of bins');
set(gca,'FontSize',16);

figure,
plot(1:nIter,timeCmp(:,3),'bo-',1:nIter,timeCmp(:,4),'rs-',1:nIter,timeCmp(:,5),'k*--');legend('ILP','IQP','Initial');
hold on;
title(sprintf('n = %d, k = %d', n,k));
set(gca,'FontSize',16);
xlabel('sample #');ylabel('function value');



