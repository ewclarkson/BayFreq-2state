%% Compare the Bayesian and frequentist approaches

%
% numExp = 1; % number of data sets to analyse
nTraj = 100;         % number of trajectories
N = 3E2;            % number of displacements
tau = 5E-3;
D1 = 10;
D2 = 1; %linspace(1,4,numExp);
sigmaB1 = sqrt(2*D1*tau);
sigmaB2 = sqrt(2*D2*tau);
p12 = 0.05;          
p21 = 0.05;
fTrue = p21/(p12+p21);
data = cell(1,nTraj);

% generate a data set
for idxTrack = 1:nTraj  
    stateVec = utilF.twoState_Markov(p12, p21, N); % random state sequence
    data{idxTrack} = utilF.brownian_displacements_2d(stateVec,sigmaB1,sigmaB2);
end

% compare frequentist method with and without flipping
[estInitial,~] = utilF.Rayleigh_2state(D1, D2, fTrue, tau, data); % NOTE: how to pick initial guess?
freq_fEst = estInitial(1); freq_D1Est = estInitial(2); freq_D2Est = estInitial(3); % estimates of f,D1,D2
[p12Est, p21Est] = utilF.method_of_moments_multi(data, freq_fEst, sqrt(2*freq_D1Est*tau), sqrt(2*freq_D2Est*tau));
[p12EstNoFlip, p21EstNoFlip] = utilF.method_of_moments_noflip(data, freq_fEst, sqrt(2*freq_D1Est*tau), sqrt(2*freq_D2Est*tau));


p12Est
p12EstNoFlip




