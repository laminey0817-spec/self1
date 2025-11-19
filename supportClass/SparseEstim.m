classdef SparseEstim < EstimIn
    % SparseEstim:  Scalar estimator class with sparsity
    
    % X = X1 with prob p1
    %   = x0 with prob 1-p1
    % 由于无稀疏性，故此部分实际不起作用
    
    properties
        x0 = 0;             % location of Dirac delta
        p1 = 0.5;           % Sparsity rate
        p0 = 0.5;           % Prior Sparsity rate
        estim1;             % Base estimator when U=1
    end
    
    methods
        % Constructor
        function obj = SparseEstim(estim1, p1)
            obj.p1 = p1;
            obj.p0 = p1;
            obj.estim1 = estim1;
        end
        
        % Set method for estim1
        function set.estim1(obj, Estim1)
            % Check to ensure input Estim1 is a valid EstimIn class
            obj.estim1 = Estim1;
        end
        
        % Compute posterior mean and variance from Gaussian estimate
        function [xhat, xvar, klDivNeg] = estim(obj, rhat, rvar)
            
            % Compute the activity probabilities
            % Get log-likelihood of rhat for U=1 and U=0
            loglike1 = obj.estim1.loglikey(rhat, rvar);
            
            rvar(rvar < eps) = eps;     % for numerical stability
            loglike0 = -( log(pi) + log(rvar) + ...
                (abs(rhat - obj.x0).^2)./rvar );
            
            % Convert log-domain quantities into posterior activity
            % probabilities (i.e., py1 = Pr{X=X1 | y}, py0 = Pr{X=x0 | y})
            exparg = loglike0 - loglike1 + log(1 - obj.p1) - log(obj.p1);
            maxarg = 500;
            exparg = max(min(exparg,maxarg),-maxarg); % numerical robustness
            py1 = (1 + exp(exparg)).^(-1);
            py0 = 1 - py1;
            
            % Compute mean and variance
            [xhat1, xvar1, klDivNeg1] = obj.estim1.estim(rhat, rvar);
            
            xhat = py1.*xhat1 + py0.*obj.x0;
            xvar = py1.*(abs(xhat1).^2 - abs(xhat).^2) + py1.*xvar1 ...
                + py0.*(abs(obj.x0).^2 -abs(xhat).^2);
            
            % Compute negative K-L divergence
            klDivNeg = py1.*klDivNeg1 ...
                + py1.*log(max(1e-8,obj.p1)./max(py1,1e-8)) ...
                + py0.*log(max(1e-8,(1-obj.p1))./max(py0,1e-8));
            
        end
    end
end

