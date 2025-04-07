function out = LogLikelihood2D(D1, D2, p12, p21, deltaT, trackDisps)

    % 
    % Calculates log-likelihood for a 2d trajectory using the 
    % forward algorithm in log-likelihood form
    % 
    % Input:
    % D1 = diffusion constant for state 1
    % D2 = diffusion constant for state 2
    % p12 = transition probability from state 1 to state 2
    % p21 = transition probability from state 2 to state 1
    % deltaT = sampling time
    % trackDisps = vector of trajectory displacements
    % 
    % Output:
    % out = Exact logarithmic likelihood of the trajectory, given the 
    % input parameters 
    % 
    % Dependencies:
    % logsumexp2.m
    % 
    
    if D1 <= 0 || D2 <= 0 || p12 <= 0 || p12 >= 1 || p21 <= 0 || p21 >= 1 % these inputs are not allowed (because they are unphysical)
        
        out = -Inf; % return a log-likelihood of -Inf
    
    else % continue with the algorithm

        import('utilB.logsumexp2')

        % --------------------- define constants ------------------------------

        p11 = 1-p12; % probability to stay in state 1, implicitly determined by p12
        p22 = 1-p21; % probability to stay in state 2, implicitly determined by p21

        N = length(trackDisps); % number of displacements in the track

        pi_1 = p21/(p12+p21); % the fraction of steps that the particle is in state 1
        pi_2 = p12/(p12+p21); % the fraction of steps that the particle is in state 2


        % ----------- define functions and compute log-likelihoods -----------

        loglikelihood_disp = @(disp_j, D_i, deltaT) -log(4*pi*D_i*deltaT)-disp_j.^2/(4*D_i*deltaT); % the likelihood for every individual displacement

        likeli_1 = loglikelihood_disp(trackDisps, D1, deltaT); % compute the likelihoods for the displacements of the track, given D1
        likeli_2 = loglikelihood_disp(trackDisps, D2, deltaT); % compute the likelihoods for the displacements of the track, given D2


        % -------------------- compute all alphas's --------------------------

        alpha_11 = log(pi_1) + likeli_1(1); % compute alpha1(1)
        alpha_12 = log(pi_2) + likeli_2(1); % compute alpha1(2)
        alphaArr_1 = zeros(1,N); alphaArr_1(1) = alpha_11; % initilise array for all values of alpha_1
        alphaArr_2 = zeros(1,N); alphaArr_2(1) = alpha_12; % initilise array for all values of alpha_2

        for j = 2:N % do for step 2 to step N
            
            alphaArr_1(j) = logsumexp2(alphaArr_1(j-1) + log(p11), alphaArr_2(j-1) + log(p21)) + likeli_1(j); 
                                                                                                          
            alphaArr_2(j) = logsumexp2(alphaArr_1(j-1) + log(p12), alphaArr_2(j-1) + log(p22)) + likeli_2(j); 
                                                                                                          
        end    

        out = logsumexp2(alphaArr_1(N), alphaArr_2(N));
    end
end