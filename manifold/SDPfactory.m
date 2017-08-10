function M = SDPfactory(n, k, A, b)
% Manifold of n-by-n symmetric positive semidefinite matrices X of rank k
% satisfying trace(A_i*X) = b_i for all i.
%
% function M = SDPfactory(n, k, A, b)
%
% A point X on the manifold is parameterized as YY^T where Y is a matrix of
% size nxk. As such, X is symmetric, positive semidefinite. We restrict to
% full-rank Y's, such that X has rank exactly k. The point X is numerically
% represented by Y (this is more efficient than working with X, which may
% be big). Tangent vectors are represented as matrices of the same size as
% Y, call them Ydot, so that Xdot = Y Ydot' + Ydot Y and trace(Xdot) == 0.
% The metric is the canonical Euclidean metric on Y.
%
%
% Note that this geometry formally breaks down at rank-deficient Y's.
% As an alternative, you may use the sphere manifold (it has larger
% dimension (by 1), but does not break down at rank drop.)
%
% The geometry is taken from the 2010 paper:
% M. Journee, P.-A. Absil, F. Bach and R. Sepulchre,
% "Low-Rank Optimization on the Cone of Positive Semidefinite Matrices".
% Paper link: http://www.di.ens.fr/~fbach/journee2010_sdp.pdf
% 
% 
% Please cite the Manopt paper as well as the research paper:
%     @Article{journee2010low,
%       Title   = {Low-rank optimization on the cone of positive semidefinite matrices},
%       Author  = {Journ{\'e}e, M. and Bach, F. and Absil, P.-A. and Sepulchre, R.},
%       Journal = {SIAM Journal on Optimization},
%       Year    = {2010},
%       Number  = {5},
%       Pages   = {2327--2351},
%       Volume  = {20},
%       Doi     = {10.1137/080731359}
%     }
% 
%
% See also: spherefactory elliptopefactory  spectahedronfactory 
%           symfixedrankYYfactory

% This file is part of Manopt: www.manopt.org.
% Original author: Bamdev Mishra, July 11, 2013.
% Contributors: Nicolas Boumal
% Change log:
%
%   April 30, 2017 (NB):
%       Adapted from spectahedronfactory by Nils Smit
    
    % Number of constraints
    nCon = length(A);
    
    M.name = @() sprintf('YY'' quotient manifold of %dx%d psd matrices of rank %d with trace constraints', n, k);
    
    M.dim = @() n*k - k*(k-1)/2; % - nCon
    
    % Euclidean metric on the total space
    M.inner = @(Y, eta, zeta) eta(:)'*zeta(:);
    
    M.norm = @(Y, eta) sqrt(M.inner(Y, eta, eta));
    
    M.dist = @(Y, Z) error('spectrahedronfactory.dist not implemented yet.');
    
    M.typicaldist = @() 1;
    
    M.proj = @projection;
    
    % Define projection function
    function etaproj = projection(Y, eta)
        % Projection onto the tangent space, i.e., on the tangent space of
        % trace(A_i*Y) = b_i forall i
        etaproj = eta;
        for i = 1:nCon
            AY = A{i}*Y;
            
            alpha = eta(:).'*AY(:)/(AY(:).'*AY(:));
            etaproj = etaproj - alpha*AY;
        end
        
        % Projection onto the horizontal space
        A_lyap = Y.'*Y;
        Q_lyap = -comm(Y.',eta);
        Omega = lyap(A_lyap, Q_lyap);
        etaproj = etaproj - Y*Omega;
    end
    
    M.tangent = M.proj;
    M.tangent2ambient = @(Y, eta) eta;
    
    M.retr = @retractionPROJ;
%     function Ynew = retractionGEO(Y, eta, t)
%         if nargin < 3
%             t = 1.0;
%         end
%         Ybar = Y + t*eta;
%         Ynew = Ybar;
%         for i = 1:nCon
%             AY = A{i}*Ybar;
%             A2Y = A{i}*AY;
%             c0 = Ybar(:).'*AY(:) - b{i};
%             c1 = AY(:).'*AY(:);
%             c2 = AY(:).'*A2Y(:);
%             
%             alpha = (-c1 + sqrt(c1^2 - 4*c2*c0))/(2*c2);
%             Ynew = Ynew + alpha*A{i}*Ybar;
%         end
%     end
%   The projection version of the retraction doesn't generally work 
%   (imaginary results)
    function Ynew = retractionPROJ(Y, eta, t)
        if nargin < 3
            t = 1.0;
        end
        Ybar = Y + t*eta;
        Ynew = Ybar;
        for i = 1:nCon
            AY = A{i}*Ybar;
            A2Y = A{i}*AY;
            c0 = Ybar(:).'*AY(:) - b{i};
            c1 = AY(:).'*AY(:);
            c2 = AY(:).'*A2Y(:);
            
            if c1^2 - 4*c2*c0 < 0
                error('Stepsize too big')
            end
            
            alpha = (-c1 + sqrt(c1^2 - 4*c2*c0))/(2*c2);
            Ynew = Ynew + alpha*A{i}*Ybar;
        end
    end
    
    M.egrad2rgrad = @projection;
    
    M.ehess2rhess = @ehess2rhess;
    function Hess = ehess2rhess(Y, egrad, ehess, eta)
        YtEta = Y.'*eta;
        
        A_ly_Om = Y.'*Y;
        Q_ly_Om = -(YtEta - YtEta.');
        Omega = lyap(A_ly_Om, Q_ly_Om);
        
        A_ly_DOm = A_ly_Om;
        Q_ly_DOm = -(comm(eta.',egrad) + comm(Y.',ehess) + comm(Omega.',YtEta + YtEta.'));
        DOmega = lyap(A_ly_DOm, Q_ly_DOm);
        
        Hess = ehess  - eta*Omega - Y*DOmega;
        
        for i = 1:nCon
            AY = A{i}*Y;
            Aeta = A{i}*eta;
            trYA2Y = (AY(:).'*AY(:));
            
            alpha = eta(:).'*AY(:)/trYA2Y;
            Dalpha = 1/trYA2Y*(ehess(:).'*AY(:)+egrad(:).'*Aeta(:)) - ...
                     eta(:).'*AY(:)/trYA2Y^2 * (2*Aeta(:).'*AY(:));
                 
            Hess = Hess - alpha*Aeta - Dalpha*AY;
        end
        
        % Project on the horizontal space
        Hess = M.proj(Y, Hess);
    end
    
    M.exp = @exponential;
    function Ynew = exponential(Y, eta, t)
        if nargin < 3
            t = 1.0;
        end
        
        Ynew = retraction(Y, eta, t);
        warning('manopt:spectrahedronfactory:exp', ...
            ['Exponential for fixed rank SDP ' ...
            'manifold not implenented yet. Used retraction instead.']);
    end
    
    % Notice that the hash of two equivalent points will be different...
    M.hash = @(Y) ['z' hashmd5(Y(:))];
    
    M.rand = @random;
    
    function Y = random()
        Y = randn(n, k);
        warning( 'manopt:SDPfactory:exp', ...
                ['Random for fixed rank SDP ' ...
                 'manifold not implenented yet.']);
    end
    
    M.randvec = @randomvec;
    function eta = randomvec(Y)
        eta = randn(n, k);
        eta = projection(Y, eta);
        nrm = M.norm(Y, eta);
        eta = eta / nrm;
    end
    
    M.lincomb = @matrixlincomb;
    
    M.zerovec = @(Y) zeros(n, k);
    
    M.transp = @(Y1, Y2, d) projection(Y2, d);
    
    M.vec = @(Y, u_mat) u_mat(:);
    M.mat = @(Y, u_vec) reshape(u_vec, [n, k]);
    M.vecmatareisometries = @() true;
end

function res = comm(x,y)
    xy = x*y;
    res = xy - xy.';
end
function res = anticomm(x,y)
    xy = x*y;
    res = xy + xy.';
end