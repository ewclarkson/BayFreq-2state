function [deadPoints,thetaMLE,logZ] = nestedsampling2D(nLive, StopRatio, priorLimits, deltaT, dispCell)
    % 
    % Perform the nested sampling algorithm. Compute Bayesian estimates, 
    % evidence and posterior. 
    % 
    % Input:
    % nLive = number of live points
    % StopRatio = stopping criterion as the fractional change in estimated 
    %             evidence falls below this value
    % priorLimits = array of upper and lower boundaries for uniform priors
    % deltaT = sampling time interval
    % dispCell = cell array containing trajectories of displacements
    % 
    % Output:
    % deadPoints = samled coordinates with posterior weights
    % thetaMLE = maximum likelihood estimate of parameters
    % logZ = estimated logarithm of evidence
    % 
    % Dependencies:
    % Enslikelihood2D.m
    % LogLikelihood2D.m
    % logsumexp2.m
    % drawlivepoint.m
    % 
    
    import('utilB.EnsLikelihood2D')
    import('utilB.logsumexp2')
    import('utilB.drawlivepoint2D')

    % ----------------- set up live points ---------------------------

    % Initialise variables
    logStopRatio = log(StopRatio); % logarithm of evidence ratio
    logZ = log(0); % initialise evidence
    logZratio = log(Inf); % initialise log-ratio of estimated evidence
    j = 1; % iteration index
    
    % Sample nLive points from a uniform prior
    D1_min = priorLimits(1,1); D1_max = priorLimits(1,2);
    D2_min = priorLimits(2,1); D2_max = priorLimits(2,2);
    p12_min = priorLimits(3,1); p12_max = priorLimits(3,2);
    p21_min = priorLimits(4,1); p21_max = priorLimits(4,2);

    init_D1s = D1_min + (D1_max-D1_min).*rand(nLive,1); % initial D1-coordinates for every live point
    init_D2s = D2_min + (D2_max-D2_min).*rand(nLive,1); % initial D2-coordinates for every live point
    init_p12s = p12_min + (p12_max-p12_min).*rand(nLive,1); % initial p12-coordinates for every live point
    init_p21s = p21_min + (p21_max-p21_min).*rand(nLive,1); % initial p21-coordinates for every live point
    
    pointArr = zeros(nLive, 4); % for later drawing new live points
    for i = 1:nLive % do for every live point

       y = [init_D1s(i),init_D2s(i),init_p12s(i),init_p21s(i)]; 
       
       if y(2) > y(1)
           y = [y(2),y(1),y(4),y(3)]; % enforce D1>D2
       end
    
       livePoints(i).pos = y; % assign values to points
       pointArr(i,:) = y; % fill pointArr
    end

    
    % Compute the likelihood associated with each live point
    
    for idxPoint = 1:nLive
    
        paramArr = livePoints(idxPoint).pos; % position of current point
        livePoints(idxPoint).logL = EnsLikelihood2D(paramArr, deltaT, dispCell); % compute and store likelihood
    end
    
    % ----------------------- begin algorithm ----------------------

    while logZratio > logStopRatio % do until the stop criterion is met, i.e. until Z has been calculated
    
        % Identify the worst point
        [logLworst, ind_worst]=min([livePoints.logL]);
        
        % Compute evidence Z up to the jth iteration
        if j==1 % do for the first iteration
            logWeight = -log(nLive+1);
        else % do for all other iterations
            logWeight = logWeight+log(nLive)-log(nLive+1);
        end
        logZ = logsumexp2(logZ, logLworst+logWeight);
        
        % Sample a new point      
        [posNew, logLnew] = drawlivepoint2D(pointArr, logLworst, deltaT, dispCell, priorLimits); % new point
        
        % Store dead point
        deadPoints(j).pos = livePoints(ind_worst).pos;
        deadPoints(j).logL = logLworst;
        deadPoints(j).logWeight = logWeight;

        % Replace dead point
        pointArr(ind_worst,:) = posNew; % set position of new point
        livePoints(ind_worst).pos = posNew; % set position of new point
        livePoints(ind_worst).logL = logLnew; % set log-likelihood of new point
        j = j+1; % update iteration number
        
        % Compute remaining evidence of points, at iteration j
        logLremain = livePoints(1).logL;
        for i = 2:nLive % do for every live point
            logLremain = logsumexp2(logLremain,livePoints(i).logL);
        end
        logZremain = logWeight+logLremain; % remaining log-evidence
    
        logZratio = logZremain-logZ;
    end
    
    % Begin final correction
    logZ = logsumexp2(logZ, logZremain);

    % Compute all posterior weights    
    for i = 1:length(deadPoints) % do for every dead point
        
        deadPoints(i).postWt = exp(deadPoints(i).logWeight+...
            deadPoints(i).logL-logZ);      
    end
    
    for i = 1:nLive % do for every live point left
        
        deadPoints(end+1).postWt = exp(logWeight+...
            livePoints(i).logL-logZ);
        
        deadPoints(end).pos = livePoints(i).pos;
        deadPoints(end).logWeight = logWeight;
        deadPoints(end).logL = livePoints(i).logL;
    end

    % Find MLE of model parameters theta
    [~, ind_max]=max([livePoints.logL]);
    thetaMLE = livePoints(ind_max).pos;
end

