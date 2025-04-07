%% Number of false positives script

% Find out the PDF of the number false positives (1->2 events, number of segments)
% This is before doing the spin-flip.

% Simulation parameters
N=1E5;              % Trajectory length
sigmaB1 = 0.5;      % Std dev for dx and dy when in state 1
sigmaB2 = 0.8;      % Std dev for dx and dy when in state 2 % NOTE: switch so that D1 > D2, as we do for the Bayesian analysis

nSim = 1E3; % number of "experiments"

%% Find the PDF of the number of segments numerically

stateVec = ones(1,N); % generate a long 1-state trajectory
sThresh = utilF.determine_optimal_thresh_displ(sigmaB1,sigmaB2); % find threshold
% q = sum(raylrnd(sigmaB1,1,1E7)<sThresh)/1E7; % probability of being correctly labelled. % Note: better to calculate as CDF.
q = raylcdf(sThresh,sigmaB1); % probability of being correctly labelled
n12FPvec = zeros(1,nSim);
for i = 1:nSim
    
    data = utilF.brownian_displacements_2d(stateVec,sigmaB1,sigmaB2); % generate displacements
    % Threshold
    stateVecEst = zeros(1,length(data));
    stateVecEst(data <= sThresh) = 1;
    stateVecEst(data > sThresh) = 2;
    % Count the number of 1->2 events
    diffX = diff(stateVecEst);
    n12FP = sum(diffX == 1);
    n12FPvec(i) = n12FP;
end

% % Do the same thing for a long 2-state trajectory, i.e. estimate n{21,FP}
% stateVec = 2*ones(1,N); % generate a long 1-state trajectory
% q2 = raylcdf(sThresh,sigmaB2); % probability of being correctly labelled
% n21FPvec = zeros(1,nSim);
% for i = 1:nSim
%     
%     data = utilF.brownian_displacements_2d(stateVec,sigmaB1,sigmaB2); % generate displacements
%     % Threshold
%     stateVecEst = zeros(1,length(data));
%     stateVecEst(data <= sThresh) = 1;
%     stateVecEst(data > sThresh) = 2;
%     % Count the number of 1->2 events
%     diffX = diff(stateVecEst);
%     n21FP = sum(diffX == -1);
%     n21FPvec(i) = n21FP;
% end

%% Plot the distribution and fit theory

% rho12FPvec = n12FPvec/N; % crossing rate

% histogram(n12FPvec,'Normalization','pdf');
% hold on
% pd = fitdist(n12FPvec','Binomial','NTrials',N);
% X = min(n12FPvec):max(n12FPvec);
% Y = binopdf(X,N,pd.p);
% plot(X,Y,'LineWidth',1)

meanSim = mean(n12FPvec)
meanTheory = q*(1-q)*(N-1)

% meanSim21 = mean(n21FPvec)
% meanTheory21 = q2*(1-q2)*(N-1)

meanTheory1 = q*(1-q); % single trial
varSim = var(n12FPvec);
% varTheory = meanTheory*(1-2*meanTheory)+meanTheory^2*(N-1);
% varTheory1 = 2*meanTheory1^2*(1-meanTheory1); % single trial
varTheory = 2*meanTheory1^2*(1-meanTheory1)*(N-1);



%% Compute distribution by formula

function n12theory = num12joints(n12,q,N)

    % calculate the probability of having n12 1-2 joints.

    % INPUT:
    % n12 - number of 1-2 joints in a trajectory
    % q - probability of having a 1 at any one site
    % N - length of trajectory

    % OUTPUT:
    % n12theory - probability of obtaining n12 1-2 joints in a trajectory
    %             of length N with probability q of having a '1' at a site
    %             and probability (1-q) of having a '2' instead, where each
    %             site is independent (a Bernoulli process).

    B1 = 0;
    for k = 1:N-1
    
        if n12 > k && n12 < 0 % then the binomial coefficient is zero or undefinable
            continue
        else
            B1 = B1 + nchoosek(k,k-n12)*q/(1-q)^k*(nchoosek(N-k-1,n12-1)+nchoosek(N-k-1,n12));
        end
    end
    n12theory = (1-q)^N*B1;
end


%% test fit a binomial distribution
% 
% Ntrials = 1E2;
% Ybino = binornd(Ntrials,0.5,1,1E4);
% 
% pd2 = fitdist(Ybino','Binomial','NTrials',Ntrials);
% % [phat,pci] = binofit(Ybino,Ntrials); % Distribution-specific function
% % [phat2,pci2] = mle(Ybino,'distribution','Binomial',"NTrials",Ntrials); % Generic distribution function
% histogram(Ybino,'Normalization','pdf')
% hold on
% 
% x = 0:100;
% y = binopdf(x,Ntrials,0.5);
% % bar(x,y,1)
% plot(x,y,'LineWidth',1.5)
% 

















