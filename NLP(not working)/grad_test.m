function [ res ] = grad_test( func, n, points )
%Help function: test if the gradient valid
%   input:  func   --  a function handle which can calculate value at
%                      a n-dim point 'x' (as a column vector), and return 
%                      'f(x)' and 'gradf(x)'.                                         
%           n      --  dimension of 'x'.
%           points --  points where the gradient is tested, should have the 
%                      same size as x. Should be a n-by-m matrix, meaning m
%                      n-dim points are tested. If not specified, this
%                      function automatically generate 100 random points.
%   output: res    --  1: gradf is correct(default).
%                      0: gradf is wrong.
    if nargin<2,
        error('func and n must be specified.');
    end
    if ~exist('points', 'var') || isempty(points) % if no points given, randomly genarate 100 points
        points = 10*rand(n,1000)-0;
    end
    if size(points,1)~=n,
       error('the number of rows in points must be n'); 
    end
    m = size(points,2); % get the number of points
    res = 1; % default result
    
    % Ram's magic number
    h = 1e-6;
    thres = 1e-4;
    
    % start test if gradf is correct w.r.t f
    for i = 1:m,
       x = points(:,i);
       [f,gradf] = func(x);
       for j = 1:n,
          % create unit vector ej
          ej = 0*x;
          ej(j) = 1;
          
          % caltulate f_ph and f_mh, where h is added/subtracted on the
          % j-th element of x
          f_ph = func(x+h*ej);
          f_mh = func(x-h*ej);
          
          % check if the gradf is correct
          gradf_est = (f_ph-f_mh)/2/h;
          gradf_est-gradf(j,:)'
          if sum(abs(gradf_est-gradf(j,:)')>thres)~=0,
              disp(['gradf fails at point trans([',num2str(x'),'])'])
              res = 0;
              break;
          end
          
       end
       
       
       if res == 0,
           break;
       end
       
    end


end

