%% Markov chain Monte Carlo optimisation

% Find the maximum of a (log-likelihood) function by moving through
% the parameter space with a Metropolis-like scheme

% Dependencies:
% Markov_likelihood
% EnsLikelihood2D

function [latest, averaged, D1s, D2s, p12s, p21s, s1, s2, s3, s4] = MCMC(num_MC_steps, burnInTime, scaleArr, guessTbl, deltaT, trackCell)

    import('utilB.EnsLikelihood2D')
    %paramArr = [0.15, 0.005, 0.08, 0.012]; % pick a random initial position in parameter space
    % pick a random initial position in parameter space
%     paramArr = [Unif_rand(guessTbl(1,1),guessTbl(1,2),1), Unif_rand(guessTbl(2,1),guessTbl(2,2),1), Unif_rand(guessTbl(3,1),guessTbl(3,2),1), Unif_rand(guessTbl(4,1),guessTbl(4,2),1)]; 
    D1_min = guessTbl(1,1); D1_max = guessTbl(1,2);
    D2_min = guessTbl(2,1); D2_max = guessTbl(2,2);
    p12_min = guessTbl(3,1); p12_max = guessTbl(3,2);
    p21_min = guessTbl(4,1); p21_max = guessTbl(4,2);

    init_D1 = D1_min + (D1_max-D1_min)*rand; % initial D1-coordinates for every live point
    init_D2 = D2_min + (D2_max-D2_min)*rand; % initial D2-coordinates for every live point
    init_p12 = p12_min + (p12_max-p12_min)*rand; % initial p12-coordinates for every live point
    init_p21 = p21_min + (p21_max-p21_min)*rand; % initial p21-coordinates for every live point

    paramArr = [init_D1, init_D2, init_p12, init_p21]; % initial random position
    
    l_logLike = EnsLikelihood2D([paramArr(1), paramArr(2), paramArr(3), paramArr(4)], deltaT, trackCell);
    
    n = num_MC_steps;
    l_theta = paramArr;
    
    all_D1s = zeros(1,num_MC_steps); % store all D1s
    all_D2s = zeros(1,num_MC_steps); % store all D2s
    all_p12s = zeros(1,num_MC_steps); % store all p12s
    all_p21s = zeros(1,num_MC_steps); % store all p21s
    
    all_D1s(1) = paramArr(1); all_D2s(1) = paramArr(2); all_p12s(1) = paramArr(3); all_p21s(1) = paramArr(4); % store initial parameter values
    
    acMoves = [0,0,0,0]; % store the number of accepted moves for each parameter (D1, D2, p12, p21)
    prMoves = [0,0,0,0]; % store the total number of proposed moves for each parameter (D1, D2, p12, p21)
    
    
    for i = 1:round(n/4)
       for k = 1:4 
           
           l = 4*(i-1)+k-1; % keeps count of the current iteration
           
           deltaTheta = normrnd(0,sqrt(scaleArr(k))); % propose a displacement for the kth parameter
           propTheta = l_theta; propTheta(k) = propTheta(k) + deltaTheta; % propose a new state
           
           prMoves(k) = prMoves(k) + 1; % add one count to the number of proposed moves along this axis 
           
           prop_logLike = EnsLikelihood2D([propTheta(1), propTheta(2), propTheta(3), propTheta(4)], deltaT, trackCell); % calculate log-likelihood for proposed parameter set
           
           if prop_logLike >= l_logLike % then accept the proposed state
               
               l_theta = propTheta; % assign the current parameter set to be the proposed one
               l_logLike = prop_logLike; % assign the current log-likelihood to be the proposed one
               
               acMoves(k) = acMoves(k) + 1; % add one count to the number of accepted moves along parameter axis k
               
               all_D1s(l+2) = propTheta(1); all_D2s(l+2) = propTheta(2); all_p12s(l+2) = propTheta(3); all_p21s(l+2) = propTheta(4); % store new parameter values
          
           else % make a new trial for whether or not to accept the proposed state
               
               u = rand();
               if log(u) <= prop_logLike - l_logLike % then accept the proposed state
                   
                   l_theta = propTheta; % assign the current parameter set to be the proposed one
                   l_logLike = prop_logLike; % assign the current log-likelihood to be the proposed one
                   
                   acMoves(k) = acMoves(k) + 1; % add one count to the number of accepted moves along parameter axis k
                   
                   all_D1s(l+2) = propTheta(1); all_D2s(l+2) = propTheta(2); all_p12s(l+2) = propTheta(3); all_p21s(l+2) = propTheta(4); % store new parameter values
                   
               else
                    all_D1s(l+2) = l_theta(1); all_D2s(l+2) = l_theta(2); all_p12s(l+2) = l_theta(3); all_p21s(l+2) = l_theta(4); % store "new" parameter values
               end
           end
       end
    end
    
    est_D1 = mean(all_D1s(burnInTime:end)); est_D2 = mean(all_D2s(burnInTime:end)); est_p12 = mean(all_p12s(burnInTime:end)); est_p21 = mean(all_p21s(burnInTime:end)); % calculate averages after thermalisation 
    
    latest = l_theta; % output the latest value
    averaged = [est_D1, est_D2, est_p12, est_p21]; % output the averaged values after thermalisation
    
    D1s = all_D1s; D2s = all_D2s; % output all diffusion constant estimations
    p12s = all_p12s; p21s = all_p21s; % output all diffusion constant estimations
    
    s1 = acMoves(1)/prMoves(1); s2 = acMoves(2)/prMoves(2); s3 = acMoves(3)/prMoves(3); s4 = acMoves(4)/prMoves(4); % store acceptance rates for every parameter
    
end