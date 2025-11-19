classdef CAwgnEstimOut < EstimOut
    % CAwgnEstimOut:  CAWGN scalar output estimation function
    % Corresponds to an output channel of the form
    %   y = scale*z + CN(0, wvar)
    
    properties
        y;                 % Measured output
        wvar;              % Variance
        wvar_min = 1e-20;   % Minimum allowed value of wvar
    end
    
    methods
        % Constructor
        function obj = CAwgnEstimOut(y, wvar)
            obj = obj@EstimOut;
            obj.y = y;
            obj.wvar = wvar;
        end
        
        % Set methods
        function obj = set.y(obj, y)
            obj.y = y;
        end
        
        function obj = set.wvar(obj, wvar)
            obj.wvar = wvar;
        end
        
        % AWGN estimation function
        % Provides the posterior mean and variance of variable z
        % from an observation y = scale*z + w
        % where z = CN(phat,pvar), w = CN(0,wvar)
        % 计算Z的后验的均值方差
        function [zhat, zvar] = estim(obj, phat, pvar)
            
            % Extract quantities
            y1 = obj.y;
            
            % Compute posterior mean and variance
            wvar1 = obj.wvar;
            gain = pvar./(pvar + wvar1);
            zhat = (gain).*(y1-phat) + phat;
            zvar = wvar1.*gain;
            
        end
        
        % Compute log likelihood
        % For sum-product GAMP, compute
        %   E( log p_{Y|Z}(y|z) ) with z = CN(phat, pvar)
        % For max-sum GAMP compute
        %   log p_{Y|Z}(y|z) @ z = phat
        function ll = logLike(obj,phat,pvar)
            
            % Ensure variance is small positive number
            wvar_pos = max(obj.wvar_min, obj.wvar);

            % Compute log-likelihood
            predErr = (abs(obj.y-phat).^2 + pvar)./wvar_pos;
            ll = -(predErr); %return the values without summing
        end
        
    end
    
end

