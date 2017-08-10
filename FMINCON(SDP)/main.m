% clear all,close all,clc;
profile on
n = 5 ;             % number of custumers(n)
k = 2 ;              % capacity   
v = 2*n+1;

rng('shuffle');     % random seed: shuffle
vert = rand(2*n+1,2);

% vert = [    0.0439    0.0647;
%     0.2790    0.5799;
%     0.2910    0.5897;
%     0.8554    0.0995;
%     0.0824    0.5057;
%     0.5478    0.0044;
%     0.0588    0.8508];% n = 3,k = 1

% % vert = [    0.6423    0.9573;
%     0.2213    0.6203;
%     0.8371    0.6003;
%     0.9711    0.1726;
%     0.8464    0.0903;
%     0.5060    0.2553;
%     0.2789    0.8586;
%     0.7466    0.9111;
%     0.2369    0.6996]; % n =4, k = 2;
   
if ~exist('map.mat'),  % n = 5,k = 2
   rng('shuffle');     % random seed: shuffle
   vert = rand(2*n+1,2);
else
    load('map.mat');
end
[ para,constr ] = pre_process( n,k,vert,[],1);


% TEST_GRAD_MAIN( fun,constr,v );

% X_real = [  0     0     1     0     0     0     0;
%             0     0     0     0     0     1     0;
%             0     0     0     1     0     0     0;
%             0     0     0     0     0     0     1;
%             0     1     0     0     0     0     0;
%             0     0     0     0     1     0     0;
%             1     0     0     0     0     0     0];
% X_real = [0     0     0     0     1     0     0     0     0;
%      0     0     0     0     0     0     0     0     1;
%      0     1     0     0     0     0     0     0     0;
%      0     0     0     0     0     1     0     0     0;
%      0     0     0     1     0     0     0     0     0;
%      0     0     0     0     0     0     0     1     0;
%      0     0     1     0     0     0     0     0     0;
%      0     0     0     0     0     0     1     0     0;
%      1     0     0     0     0     0     0     0     0];% heuristic
% 
% X_real = [0     1     0     0     0     0     0     0     0
%      0     0     0     0     0     1     0     0     0
%      0     0     0     0     1     0     0     0     0
%      0     0     0     1     0     0     0     0     0
%      0     0     0     0     0     0     0     1     0
%      0     0     1     0     0     0     0     0     0
%      0     0     0     0     0     0     0     0     1
%      0     0     0     0     0     0     1     0     0
%      1     0     0     0     0     0     0     0     0];%IP
% 
X_real = [     0     0     0     0     0     1     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0     1;
     0     0     1     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     1     0     0     0;
     0     0     0     1     0     0     0     0     0     0     0;
     0     1     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     1     0     0     0     0;
     0     0     0     0     0     0     0     0     1     0     0;
     0     0     0     0     1     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     1     0;
     1     0     0     0     0     0     0     0     0     0     0];% n = 5,k = 2


% suppose x=[info_for_permutation_matrix; rank-1_condition_corresponding_eigenvector]
fun = @(x)quadcost(x,para.Q,v^2);
% constraint = @(x)tsp_constr_new(x,constr,v^2);
% x0 = [reshape(X_real',v^2,1)+0.5*rand(v^2,1);1e+0*rand(v^2+1,1)]; 
% x0 = rand(2*v^2+1,1);
% lb = [zeros(v^2,1);-0.5*ones(v^2+1,1)];
% ub = [ones(v^2,1);2*ones(v^2+1,1)];


constraint = @(x)tsp_constr(x,constr);
% X_real = [zeros(v-1,1),eye(v-1);1,zeros(1,v-1)];
x0 = [reshape(X_real',v^2,1)+0.5*rand(v^2,1)];
% x0 = rand(v^2,1);
lb = zeros(v^2,1);
ub = ones(v^2,1);




options = optimoptions(@fmincon,'Algorithm','interior-point','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
% % options = optimoptions(@fmincon,'Algorithm','sqp','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
% options = optimoptions(@fmincon,'Algorithm','active-set','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
% % options = optimoptions(@fmincon,'Algorithm','trust-region-reflective','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);

options.OptimalityTolerance = 1e-5;  
options.ConstraintTolerance = 1e-4;
options.StepTolerance = 1e-6;
options.MaxFunctionEvaluations = 1e+6;
options.MaxIterations = 1e+6;
[x,fval,eflag,output,lambda,grad] = fmincon(fun,x0,...
    [],[],[],[],lb,ub,constraint,options);

% A = [];B = [];Aeq = [zeros(v*(v-1),1);1;zeros(v-1,1)]';Beq = [1];
% for i = 1:length(constr.ineq.A);
%     A = [A;constr.ineq.A{i}];
%     B = [B;constr.ineq.B{i}];
% end
% for i = 1:length(constr.eq.A);
%     Aeq = [Aeq;constr.eq.A{i}];
%     Beq = [Beq;constr.eq.B{i}];
% end
% [x,fval,eflag,output,lambda,grad] = fmincon(fun,x0,...
%     A,B,Aeq,Beq,lb,ub,constraint,options);


profile viewer
profile off