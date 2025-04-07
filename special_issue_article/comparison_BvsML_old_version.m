%% Compare the Bayesian and maximum likelihood methods

% simulation parameters

%% generate and analyse data

numExp = 2; % number of data sets to analyse

nTraj = 10;         % number of trajectories
N = 3E2;            % number of sampling time intervals
tau = 5E-3;
D1vec = linspace(10,15,numExp);
D2vec = linspace(1,1,numExp);
% sigmaB1vec = sqrt(2*D1vec*tau);
% sigmaB2vec = sqrt(2*D2vec*tau);
p12vec = linspace(0.1,0.1,numExp);          
p21vec = linspace(0.1,0.1,numExp);     
dispCell = cell(1,nTraj);
bayesCell = cell(1,nTraj); % hold parameter estimates
freqCell = cell(1,nTraj);

for idx = 1:numExp

    D1 = D1vec(idx); D2 = D2vec(idx);
    p12 = p12vec(idx); p21 = p21vec(idx);
    sigmaB1 = sqrt(2*D1*tau); sigmaB2 = sqrt(2*D2*tau);

    % generate a data set
    for idxTrack = 1:nTraj  
        stateVec = utilF.twoState_Markov(p12, p21, N); % random state sequence
        dispCell{idxTrack} = utilF.brownian_displacements_2d(stateVec,sigmaB1,sigmaB2);
    end
    
    % parameters for Bayesian analysis
    % priors
    D1_max = D1+0.5*D1; % known upper limit of D1
    D1_min = max(D1-0.5*D1,0.01); % known lower limit of D1
    D2_max = D2+0.5*D2; % known upper limit of D2
    D2_min = max(D2-0.5*D2,0.01); % known lower limit of D2
    p12_max = min(5*p12,1); % known upper limit of p12
    p12_min = min(1/5*p12,0.01); % known lower limit of p12
    p21_max = min(5*p21,1); % known upper limit of p21
    p21_min = min(1/5*p21,0.01); % known lower limit of p21 
    % algorithm parameters
    nLive = 100; % number of live points
    StopRatio = 1E-4; % stop criterion for evidence
    priorLimits = [D1_min,D1_max,D2_min,D2_max,p12_min,p12_max,p21_min,p21_max];

    [finalSeq, thetaMLE, logZ] = utilB.nestedsampling2D(nLive, StopRatio, priorLimits, deltaT, dispCell);
    [logZ_error, thetaBayes, errorBayes] = utilB.bayesianestimate(finalSeq,logZ,nLive,false);
    bayesCell{idx} = thetaBayes;

    % do ML-analysis

end

%% parameters for ML analysis








%% Plot and compare estimates

% D1,D2 (=10,1 in placeholder code below)
xGT = linspace(D2,4*D2,4);
yGT = D1*ones(1,4);
plot(xGT,yGT,'o','MarkerFaceColor','#D95319','MarkerEdgeColor','none')
hold on

x = D2/3+linspace(D2,4*D2,4);
y = D1*ones(1,4);
yneg = 1E-1*[1 3 4 3];
ypos = yneg; % = [2 5 3 5 2 5 2 2 5 5];
xneg = 1E-2*[2 5 3 4];
xpos = xneg; % [2 5 3 5 2 5 2 2 5 5];
errorbar(x,y,yneg,ypos,xneg,xpos,'.','Color','#0072BD')
hold on
errorbar(x,rand(1,4).*y,yneg,ypos,xneg,xpos,'.','Color','#77AC30')
grid on
xlim([0.1,5])
ylim([0,2*D1])
xlabel('D_2')
ylabel('D_1')



















