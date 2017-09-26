function [ ] = simulation( nlist,elist,k )
% this file is used to simulate our dateset with n = nlist where n is the 
% number of customers. elist represent the index of cases that we want to 
% simulate for each n in nlist. k is the capacity with 2 as the default. 
    if nargin==2,
        k = 2; 
    elseif nargin<2,
        error('Check Input!');
    end
    
    load '../dataset.mat';
    
    len_n = length(nlist); % length of nlist
    len_e = length(elist); % length of elist
    
    %%%%%%%%%%%%% generate solver options for YALMIP %%%%%%%%%%%%%%%%
    ops = sdpsettings('verbose',2,'debug',1,'savesolveroutput',1);
    relative_gap = 0.035;
    maxtime = 7200;

    % setting cplex
    ops.cplex.mip.tolerances.mipgap=relative_gap;
    ops.cplex.MaxTime = maxtime;

    % setting gurobi
    ops.gurobi.MIPGap = relative_gap;
    ops.gurobi.TimeLimit = maxtime;

    % setting mosek
    ops.mosek.MSK_DPAR_MIO_MAX_TIME = maxtime;
    ops.mosek.MSK_DPAR_MIO_TOL_REL_GAP = relative_gap;



    % setting baron
%     ops.baron.cplexlibname = 'cplex1271.dll';  % probably we should not use cplex here
    ops.baron.epsr = relative_gap;
    ops.baron.maxtime = maxtime;

    % setting xpress
    ops.xpress.MIPRELSTOP = relative_gap;
    ops.xpress.MAXTIME = maxtime;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:len_n,
       n = nlist(i);
       v = 2*n+1;
       for j = 1:len_e,
           sample_num = elist(j);
           vert = pos{n}{sample_num}; % get map
           
           % generate cost matrix and make it psd
           c = generate_C( vert,v );
           
           % solve IQP 
           generateIQP( c,v,k,ops,n,sample_num);
           
           % solve ILP
           if n<=12, % experiment has shown that when n>12, no solver can generate feasible result for ILP in 2 hours
                generateILP( c,v,k,ops,n,sample_num);
           end
           
       end                
    end
    
%     junk = 1
end

