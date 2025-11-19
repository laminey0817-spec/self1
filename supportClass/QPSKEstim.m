classdef QPSKEstim < EstimIn
    % QPKSEstimIn:  QPSK Signal input estimation function
    % Construct: QPSKEstimIn(prob1,probm1,probi,probmi)
    
    properties
        %Probability of 1,-1,1i,-1i
        prob1 = 0.25;           
        probm1 = 0.25;
        probi = 0.25;
        probmi = 0.25;
        mean0 = 0;         % Prior mean
        var0 = 1;          % Prior variance
    end
   
 
    methods
        % Constructor
        function obj = QPSKEstim()
            obj = obj@EstimIn;
            obj.prob1 = 0.25;           
            obj.probm1 = 0.25;
           obj.probi = 0.25;
           obj.probmi = 0.25; 
            obj.mean0 = ((obj.prob1 - obj.probm1) + 1i * (obj.probi - obj.probmi))/4;
            obj.var0 = obj.prob1 .* abs( 1 - obj.mean0).^2 ...
                        + obj.probm1 .* abs( -1 - obj.mean0).^2 ...
                        + obj.probi .* abs( 1i - obj.mean0).^2 ...
                         + obj.probmi .* abs( -1i - obj.mean0).^2;
        end

        % Prior mean and variance
        function [mean0, var0, valInit] = estimInit(obj)
            mean0 = obj.mean0;
            var0  = obj.var0;
            valInit = 0;
        end
       
        % Size
        function [nx,ncol] = size(obj)
            [nx,ncol] = size(obj.mean0);
        end

        % QPSK estimation function
        % Provides the mean and variance of a variable x = CN(uhat0,uvar0)
        % from an observation rhat = x + w, w = CN(0,rvar)
        function [xhat, xvar, val] = estim(obj, rhat, rvar)
            
            % Get prior
            global alphabet;
            Xbet = length(alphabet);
            prior_prob = 1/Xbet * ones(size(rhat,1), size(rhat,2), Xbet);
            con_prob = zeros(size(rhat,1), size(rhat,2), Xbet);
            post_prob = zeros(size(rhat,1), size(rhat,2), Xbet);
            
            for Xi = 1:Xbet
                con_prob(:,:,Xi) = exp(-abs(alphabet(Xi) - rhat).^2 ./ rvar) / pi ./ rvar;
            end
            
            % 归一化
            sum_con = sum(con_prob,3); 
            sum_con((sum_con < 10^(-50))) = 10^(-50); 
            for Xi = 1:Xbet
                con_prob(:,:,Xi) = con_prob(:,:,Xi) ./ sum_con;
            end
            
            % Compute posterior mean and variance
            for Xi = 1:Xbet
                post_prob(:,:,Xi) = prior_prob(:,:,Xi) .* con_prob(:,:,Xi);
            end
            
            % 归一化
            sum_post = sum(post_prob,3);
            sum_post((sum_post < 10^(-50))) = 10^(-50);
            for Xi = 1:Xbet
                post_prob(:,:,Xi) = post_prob(:,:,Xi) ./ sum_post;
            end
            
            % Compute Xhat and Xvar
            xhat = zeros(size(rhat));
            for Xi = 1:Xbet
                xhat = xhat + post_prob(:,:,Xi) * alphabet(Xi);
            end
            xvar = max(1e-4, 1-abs(xhat).^2); % E(X^2)=1
       
            % KL散度
            uhat0 = obj.mean0; %最原始的均值和方差，是用来算散度的
            uvar0 = obj.var0;
             
            % Compute the negative KL divergence
            %   klDivNeg = \sum_i \int p(x|r)*\log( p(x) / p(x|r) )dx
            xvar_over_uvar0 = rvar./(uvar0+rvar);
            val =  (log(xvar_over_uvar0) + (1-xvar_over_uvar0) ...
                - abs(xhat-uhat0).^2./uvar0 );

         end
        
        % Computes the log-likelihood, log p(rhat), for rhat = x + v, where 
        % x = CN(obj.mean0, obj.var0) and v = CN(0,rvar)
        function logpy = loglikey(obj, rhat, rvar)
            logpy = -( log(pi) + log(obj.var0 + rvar) + ...
                (abs(rhat - obj.mean0).^2) ./ (obj.var0 + rvar) );
        end    
    end
    
end

