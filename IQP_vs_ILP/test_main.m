clear all;close all;clc

nIter = 100;
n = 8;
k = 1;
appx = 1e-04;
for i = 1:nIter
    i
    [outputILP,obj_ILP,solution_ILP,outputIQP,obj_IQP,solution_IQP] = main_hp_func(n,k,i,appx);
    timeCmp(i,1) = outputILP.solvertime;
    timeCmp(i,2) = outputIQP.solvertime;
    timeCmp(i,3) = obj_ILP;
    timeCmp(i,4) = obj_IQP;
end


figure,
plot(1:nIter,timeCmp(:,1),'bo-',1:nIter,timeCmp(:,2),'rs-');legend('ILP','IQP');
hold on;
title(sprintf('n = %d, k = %d', n,k));
set(gca,'FontSize',16);
xlabel('sample #');ylabel('solver time');
print -dpng


figure,
hist(timeCmp(:,1:2));hold on;
title(sprintf('n = %d, k = %d', n,k));
legend('ILP','IQP');
xlabel('solver time');
ylabel('number of bins');
set(gca,'FontSize',16);
print -dpng

figure,
plot(1:nIter,timeCmp(:,3),'bo-',1:nIter,timeCmp(:,4),'rs-');legend('ILP','IQP');
hold on;
title(sprintf('n = %d, k = %d', n,k));
set(gca,'FontSize',16);
xlabel('sample #');ylabel('solver time');
print -dpng
