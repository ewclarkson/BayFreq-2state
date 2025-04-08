function out = EnsLikelihood2D(paramArr, deltaT, dispCell)

    % 
    % Calculates the log-likelihood for an ensemble of 2d trajectories.
    % Trajectories are assumed to be independent.
    % 
    % Input:
    % paramArr = a position in parameter space [D1, D2, p12, p21]
    % deltaT = sampling time (time between frames)
    % dispCell = cell containing all individual displacements, on the form 
    %            {[s1,s2,...], [r1,r2,...],...}
    % 
    % Output:
    % logLtot = total log-likelihood of all trajectories
    % 
    % Dependencies:
    % LogLikelihood2D
    % 

    logLtot = 0; % initialise accumulated log-likelihood
    for n = 1:length(dispCell) % do for every sequence
        
       loglikeli_n = utilB.LogLikelihood2D(paramArr(1), paramArr(2), paramArr(3), paramArr(4), deltaT, dispCell{n});
       logLtot = logLtot + loglikeli_n; % add this log-likelihood to the sum of all such
    end    

    out = logLtot; % function output
end