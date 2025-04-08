function fVal = func_self_consistency_eq( p , f, D12, lThresh1X , lThresh2X)
    %
    % Function whose root is an estimate of p_12.
    %
    % Input:
    % p = argument
    % f = fraction of time spent in state 1
    % D12 = (no of 1-2 transitions - n12FalsePositivesEst)/
    %        (time spent in state 1)
    % lThresh1X = length threshold for segments of type 1X
    % lThresh2X = length threshold for segments of type 2X 
    %
    % Output:
    % fVal = value of function
    %

    fVal = D12 + p.*(1-(1-p).^lThresh1X) + p.*(1-(1-(f/(1-f))*p).^lThresh2X) - p; % assuming long trajectories
    % fVal = D12 + p.*( (1-(1-p).^lThresh1X)./(1-(1-p).^N) + (1-(1-(f/(1-f))*p).^lThresh2X)./(1-(1-(f/(1-f))*p).^N) ) - p; % with proper normalisation

end

