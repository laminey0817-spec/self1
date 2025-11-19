classdef BiGAMPOpt
    % Options for the BiG-AMP optimizer.
    
    properties

        %Return additional outputs useful for EM learning. Can be disabled
        %to reduce the memory required for the output structure
        saveEM = true;

        %Number of iterations
        nit = 1500;
        
        %Minimum number of iterations- sometimes useful with warm starting
        %approaches
        nitMin = 30; %0 for no effect
        
        %Specify a convergence tolerance. Iterations are terminated when
        %the norm of the differnece in two iterates divided by the norm of
        %the current iterate falls below this threshold. Set to -1 to
        %disable
        tol = 1e-8;
                
        %***** Initialization
        xhat0 = [];
        Ahat0 = [];
        xvar0 = [];
        Avar0 = [];

        pvarStep = true;
        
        %Initial step size, or fixed size for non-adaptive steps
        step = 0.05;
        
        % Minimum step size
        stepMin = 0.05;
        
        % Maximum step size
        stepMax = 0.5;
        
        % Multiplicative step size increase, when successful
        stepIncr = 1.1;
        
        % Multiplicative step size decrease, when unsuccessful
        stepDecr = 0.5;
        
        %Maximum number of allowed failed steps before we decrease stepMax,
        %inf for no effect
        maxBadSteps = inf;
        
        %Amount to decrease stepMax after maxBadSteps failed steps, 1 for
        %no effect
        maxStepDecr = 0.5;
        
        % Iterations are termined when the step size becomes smaller
        % than this value. Set to -1 to disable
        stepTol = -1;

        %***** Variance Control
        %Minimum variances. See code for details of use.
        pvarMin = 1e-13;
        xvarMin = 0;
        AvarMin = 0;
        zvarToPvarMax = 0.99;   % prevents negative svar, should be near 1
        
        %Variance threshold for rvar and qvar, set large for no effect
        varThresh = 1e6;
    end
    
    methods
        
        % Constructor with default options
        function opt = BiGAMPOpt(varargin)
            return
        end
    end
    
end
