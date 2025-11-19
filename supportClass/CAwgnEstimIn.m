classdef CAwgnEstimIn < EstimIn
    % CAwgnEstimIn:  Circular AWGN scalar input estimation function
    
    properties
        var0_min = eps;    % Minimum allowed value of var0
        mean0 = 0;         % Prior mean
        var0 = 1;          % Prior variance
    end
   
    properties (Hidden)
        mixWeight = 1;              % Weights for autoTuning
    end
 
    methods
        % Constructor
        function obj = CAwgnEstimIn(mean0, var0)
            obj = obj@EstimIn;
            obj.mean0 = mean0;
            obj.var0 = var0;
        end

        %Set Methods     
        function obj = set.mean0(obj, mean0)
            obj.mean0 = mean0;
        end
        
        function obj = set.var0(obj, var0)
            obj.var0 = max(obj.var0_min,var0); % avoid too-small variances!
        end

        % Prior mean and variance
        function [mean0, var0, valInit] = estimInit(obj)
            mean0 = obj.mean0;
            var0  = obj.var0;
            valInit = 0;
        end

        % Circular AWGN estimation function
        % Provides the mean and variance of a variable x = CN(uhat0,uvar0)
        % from an observation rhat = x + w, w = CN(0,rvar)
        function [Ahat, Avar, val] = estim(obj, qhat, qvar)
            % Get prior
            uhat0 = obj.mean0;
            uvar0 = obj.var0; 
            
            % Compute posterior mean and variance
            gain = uvar0./(uvar0+qvar);
            Ahat = gain.*(qhat-uhat0)+uhat0;
            Avar = gain.*qvar;
            
            % Compute the negative KL divergence
            %   klDivNeg = \sum_i \int p(x|r)*\log( p(x) / p(x|r) )dx
            xvar_over_uvar0 = qvar./(uvar0+qvar);
            val =  (log(xvar_over_uvar0) + (1-xvar_over_uvar0) ...
                - abs(Ahat-uhat0).^2./uvar0 );

        end
        
        % Computes the likelihood p(rhat) for rhat = x + v, v = CN(0,rvar)
        function py = plikey(obj,rhat,rvar)
            py = exp(-1./((obj.var0+rvar)).*abs(rhat-obj.mean0).^2);
            py = py./ (pi*(obj.var0+rvar));
        end
        
        % Computes the log-likelihood, log p(rhat), for rhat = x + v, where 
        % x = CN(obj.mean0, obj.var0) and v = CN(0,rvar)
        function logpy = loglikey(obj, rhat, rvar)
            logpy = -( log(pi) + log(obj.var0 + rvar) + ...
                (abs(rhat - obj.mean0).^2) ./ (obj.var0 + rvar) );
        end
        
    end
    
end

