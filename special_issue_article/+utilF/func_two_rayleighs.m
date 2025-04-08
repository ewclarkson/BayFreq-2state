function fVal = func_two_rayleighs(s, sigmaMin, sigmaMax)
    %
    % Function which depends on two Rayleigh CDFs. 
    %
    % Input:
    % s = argument
    % sigmaMin = smallest scale parameter
    % sigmaMax = largest scale parameter
    %
    % Output:
    % fVal = value of function
    %

       fVal = 1 - raylcdf(s,sigmaMin) - raylcdf(s,sigmaMax);  
  
end

