function estFin = BiGAMP(gX, gA, gOut, opt)
% BiGAMP:  Bilinear Generalized Approximate Message Passing
%   X,A -> Z = A*X -> Y,
% where the components of X and A are independent and the mapping Z -> Y is
% separable. X is NxL, A is MxN, and Z,Y are consequently MxL.

% INPUTS:
% gX:  An input estimator derived from the EstimIn class
%    based on the input distribution p_x_{nl}(x_nl).
% gA:  An input estimator derived from the EstimIn class
%    based on the input distribution p_a_{mn}(a_mn).
% gOut:  An output estimator derived from the EstimOut class
%    based on the output distribution p_{Y|Z}(y_ml|z_ml).
% opt (optional):  A set of options of the class BiGAMPOpt.
% OUTPUTS:
% estFin: Structure containing final BiG-AMP outputs

%% Setup
nit     = opt.nit;              % number of iterations
nitMin  = opt.nitMin;           % minimum number of iterations
step    = opt.step;             % step size
stepMin = opt.stepMin;          % minimum step size
stepMax = opt.stepMax;          % maximum step size
stepIncr = opt.stepIncr;        % step inc on succesful step
stepDecr = opt.stepDecr;        % step dec on failed step
tol = opt.tol;                  % Convergence tolerance
stepTol = opt.stepTol;          % minimum allowed step size
pvarStep = opt.pvarStep;        % incldue step size in pvar/zvar
maxBadSteps = opt.maxBadSteps;  % maximum number of allowed bad steps
maxStepDecr = opt.maxStepDecr;  % amount to decrease maxStep after failures
zvarToPvarMax = opt.zvarToPvarMax;  % maximum zvar/pvar ratio

%Determine requested outputs
saveEM = opt.saveEM;

%Assign Avar and xvar mins
xvarMin = opt.xvarMin;
Avarmin = opt.AvarMin;

%% Initialization

%Replace these defaults with the warm start values if provided in the
%options object
if ~isempty(opt.xhat0)
    xhat = opt.xhat0;
    valIn = -inf;
end
if ~isempty(opt.xvar0)
    xvar = opt.xvar0;
end
if ~isempty(opt.Ahat0)
    Ahat = opt.Ahat0;
end
if ~isempty(opt.Avar0)
    Avar = opt.Avar0;
end

%Placeholder initializations- values are not used
xhatBar = 0;
AhatBar = 0;
shat = 0;
svar = 0;
pvarOpt = 0;
zvarOpt = 0;

%Init valOpt empty
valOpt = [];

%Specify minimum variances
pvarMin = opt.pvarMin;

%Placeholder initializations
rhat = 0;
rvar = 0;
qhat = 0;
qvar = 0;

%Cost init
val = zeros(nit,1);
zhatOpt = 0;

%% Iterations
%Control variable to end the iterations
stop = false;
it = 0;
failCount = 0;

%Handle first step
step1 = 1;

% Main iteration loop
while ~stop
    
    % Iteration count
    it = it + 1;
    
    % Check for final iteration
    if it >= nit
        stop = true;
    end
    
    %（R1）
    Ahat2 = abs(Ahat).^2;
    xhat2 = abs(xhat).^2;
    
    %Compute zvar (R1)
    zvar = Avar*xhat2 + Ahat2*xvar;
    
    % (R3)
    pvar = zvar + Avar*xvar;
    
    %Include pvar step
    if pvarStep
        pvar = step1*pvar + (1-step1)*pvarOpt;
        zvar = step1*zvar + (1-step1)*zvarOpt;
    end
    
    %Update zhat (R2)
    zhat = Ahat * xhat;
    
    % Continued output step (R4)
    phat = zhat- shat.*zvar;
    
    % Compute log likelihood at the output and add it the total negative
    % K-L distance at the input.
    valOut = sum(sum(gOut.logLike(zhat,pvar)));
    val(it) = valOut + valIn;
    
    % Determine if candidate passed
    if ~isempty(valOpt)
        
        stopInd = length(valOpt);
        startInd = max(1,stopInd - 1);
        
        %Check the step
        pass = (val(it) > min(valOpt(startInd:stopInd))) || (step <= stepMin);
        
    else
        pass = true;
    end

    % If pass, set the optimal values and compute a new target shat and
    % snew.
    if (pass)
        
        %Slightly inrease step size after pass
        step = stepIncr*step;
        
        % Set new optimal values
        shatOpt = shat;
        svarOpt = svar;
        xhatBarOpt = xhatBar;
        xhatOpt = xhat;
        AhatBarOpt = AhatBar;
        AhatOpt = Ahat;
        pvarOpt = pvar;
        zvarOpt = zvar;
        
        %Bound pvar
        pvar = max(pvar, pvarMin);
        
        %We keep a record of only the succesful step valOpt values
        valOpt = [valOpt val(it)]; %#ok<AGROW>
        
        % Output nonlinear step
        % (R5)(R6)
        [zhat0,zvar0] = gOut.estim(phat,pvar);
            
        %Compute 1/pvar
        pvarInv = 1 ./ pvar;
        
        %Update the shat quantities
        %（R8）步骤
        shatNew = pvarInv.*(zhat0-phat);
        %（R7）步骤 min是for数据的稳定性
        svarNew = pvarInv.*(1-min(zvar0./pvar,zvarToPvarMax));
   
        %Enforce step size bounds
        step = min([max([step stepMin]) stepMax]);
        
    else
        
        %Check on failure count
        failCount = failCount + 1;
        if failCount > maxBadSteps
            failCount = 0;
            stepMax = max(stepMin,maxStepDecr*stepMax);
        end
        % Decrease step size
        step = max(stepMin, stepDecr*step);
        
        %Check for minimum step size
        if step < stepTol
            stop = true;
        end
    end
    
    % Check for convergence if step was succesful
    if pass
        testVal = norm(zhat(:) - zhatOpt(:)) / norm(zhat(:));
        if (it > 1) && ...
                (testVal < tol)
            stop = true;
        end
        % 循环终止条件
        
        %Set other optimal values- not actually used by iterations
        AvarOpt = Avar;
        xvarOpt = xvar;
        zhatOpt = zhat;
        
        %Save EM variables if requested
        if saveEM
            rhatFinal = rhat;
            rvarFinal = rvar;
            qhatFinal = qhat;
            qvarFinal = qvar;
            zvarFinal = zvar0;
            pvarFinal = pvar;
        end
    end
    
    % Create new candidate shat
    if it > 1 
        step1 = step;
    end
    
    shat = (1-step1)*shatOpt + step1*shatNew;
    svar = (1-step1)*svarOpt + step1*svarNew;
    xhatBar = (1-step1)*xhatBarOpt + step1*xhatOpt;
    AhatBar = (1-step1)*AhatBarOpt + step1*AhatOpt;
    
    %Compute rvar and correct for infinite variance
    % （R9）步骤 计算X的似然分布的方差
    rvar = 1./((abs(AhatBar).^2)'*svar);
    rvar(rvar > opt.varThresh) = opt.varThresh;
    
    %Update rhat
    % （R10）步骤中的一部分
    rGain = (1 - (rvar.*(Avar'*svar)));
    
    rGain = min(1,max(0,rGain));
    % 计算X的似然分布的均值
    rhat = xhatBar.*rGain +...
        rvar.*(AhatBar'*shat);
    rvar = max(rvar, xvarMin);
    
    % （R11）计算A的似然分布的方差
    qvar = 1./(svar*(abs(xhatBar).^2)');
    qvar(qvar > opt.varThresh) = opt.varThresh;
    
    % Update qhat
    % （R12）
    qGain = (1 - (qvar.*(svar*xvar')));
    qGain = min(1,max(0,qGain));
    qhat = AhatBar.*qGain +...
        qvar.*(shat*xhatBar');
    qvar = max(qvar,Avarmin);
    
    % Input nonlinear step
    % (R13) - (R16)
    [Ahat,Avar,valInA] = gA.estim(qhat, qvar); % 得到信道的后验
    [xhat,xvar,valInX] = gX.estim(rhat, rvar); % 得到信号的后验
    
    %Update valIn
    valIn = sum( valInX(:) ) + sum ( valInA(:) );
    
    %Don't stop before minimum iteration count
    if it < nitMin
        stop = false;
    end
    
end

%% Save the final values

%Estimates of the two matrix factors
estFin.xhat = xhatOpt;
estFin.xvar = xvarOpt;
estFin.Ahat = AhatOpt;
estFin.Avar = AvarOpt;
estFin.rhat = rhatFinal;
estFin.rvar = rvarFinal;
estFin.qhat = qhatFinal;
estFin.qvar = qvarFinal;
estFin.zhat = zhatOpt;
estFin.zvar = zvarFinal;
estFin.phat = phat;
estFin.pvar = pvarFinal;
estFin.BiGcycletimes = it;

