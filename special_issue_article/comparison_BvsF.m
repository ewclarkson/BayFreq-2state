%% Compare the Bayesian and frequentist approaches

% Analyse simulated data with ns-HMM and tpe-RMM.
%
% ns-HMM: a two-state hidden Markov model with nested sampling for 
%         sampling posterior distribution
% tpe-RMM: a Rayleigh mixture model with added estimation of
%          transition probabilities

% ----------- Set data parameters -----------------
numExp = 4; % number of data sets to analyse
nTraj = 100;         % number of trajectories
N = 3E2;            % number of displacements
tau = 5E-3;
D1vec = [10,10,10,10]; %linspace(10,10,numExp);
D2vec = [1,3,1,3]; %linspace(1,4,numExp);
p12vec = [0.01 0.01 0.05 0.05]; %linspace(0.01,0.01,numExp);          
p21vec = [0.01 0.01 0.05 0.05]; %linspace(0.01,0.01,numExp);
% fTrueVec = p21vec./(p12vec+p21vec);
data = cell(1,nTraj);
% -------------------------------------------------

BayesEst = cell(1,numExp); % hold bayesian parameter estimates
BayesError = cell(1,numExp); % hold bayesian error estimates
BayesSeq = cell(1,numExp); % full NS output
freqEst = cell(1,numExp); % hold frequentist parameter estimates

for idx = 1:numExp
   
    D1 = D1vec(idx); D2 = D2vec(idx);
    p12 = p12vec(idx); p21 = p21vec(idx); fTrue = p21/(p12+p21);
    sigmaB1 = sqrt(2*D1*tau); sigmaB2 = sqrt(2*D2*tau);

    % generate a data set
    for idxTrack = 1:nTraj  
        stateVec = utilF.twoState_Markov(p12, p21, N); % random state sequence
        data{idxTrack} = utilF.brownian_displacements_2d(stateVec,sigmaB1,sigmaB2);
    end

    % ----------- do bayesian analysis ---------------
    D1_max = 20; %D1+0.5*D1; % known upper limit of D1
    D1_min = 3; %max(D1-0.5*D1,0.01); % known lower limit of D1
    D2_max = 8; %D2+0.5*D2; % known upper limit of D2
    D2_min = 0.3; %max(D2-0.5*D2,0.01); % known lower limit of D2
    p12_max = 0.2; %min(5*p12,1); % known upper limit of p12
    p12_min = 0.005; %min(1/5*p12,0.01); % known lower limit of p12
    p21_max = 0.2; %min(5*p21,1); % known upper limit of p21
    p21_min = 0.005; %min(1/5*p21,0.01); % known lower limit of p21 
    % algorithm parameters
    nLive = 300; % number of live points
    StopRatio = 1E-4; % stop criterion for evidence
    priorLimits = [D1_min,D1_max;D2_min,D2_max;p12_min,p12_max;p21_min,p21_max];

    tic
    [finalSeq, thetaMLE, logZ] = utilB.nestedsampling2D(nLive, StopRatio, priorLimits, tau, data);
    [logZ_error, thetaBayes, errorBayes] = utilB.bayesianestimate(finalSeq,logZ,nLive,false);
    BayesEst{idx} = thetaBayes; % store parameter estimates
    BayesError{idx} = errorBayes; % store error estimates
    BayesSeq{idx} = finalSeq; % store the full sequence
    toc

    % ----------- do frequentist analysis ---------------
    tic
    [estInitial,~] = utilF.Rayleigh_2state((D1_min+D1_max)/2, (D2_min+D2_max)/2, (p21_max/(p21_max+p12_min)+p21_min/(p21_min+p12_max))/2, tau, data); 
    freq_fEst = estInitial(1); freq_D1Est = estInitial(2); freq_D2Est = estInitial(3); % estimates of f,D1,D2
    [freq_p12Est, freq_p21Est] = utilF.wrapperFreq(data, freq_fEst, sqrt(2*freq_D1Est*tau), sqrt(2*freq_D2Est*tau));
    freqEst{idx} = [freq_D1Est,freq_D2Est,freq_p12Est,freq_p21Est]; % store parameter estimates
    toc
end

% ------------- extract and plot estimates ---------------

BayesPos = cell(1,numExp); % holds all positions
BayesPW = cell(1,numExp); % holds all posterior weights

D1freqEst = zeros(1,numExp); %D1MLerror = zeros(1,numExp);
D2freqEst = zeros(1,numExp); %D2MLerror = zeros(1,numExp);
p12freqEst = zeros(1,numExp); %p12MLerror = zeros(1,numExp);
p21freqEst = zeros(1,numExp); %p21MLerror = zeros(1,numExp);

for idx = 1:numExp
    BayesPosTemp = [];
    BayesPWtemp = [];
    for idy = 1:length(BayesSeq{idx})
        BayesPosTemp(idy,:) = BayesSeq{idx}(idy).pos;
        BayesPWtemp(idy) = BayesSeq{idx}(idy).postWt;
    end
    BayesPos{idx} = BayesPosTemp;
    BayesPW{idx} = BayesPWtemp;

    D1freqEst(idx) = freqEst{idx}(1); %D1MLerror(idx) = MLerror{idx}(1);
    D2freqEst(idx) = freqEst{idx}(2); %D2MLerror(idx) = MLerror{idx}(2);
    p12freqEst(idx) = freqEst{idx}(3); %p12MLerror(idx) = MLerror{idx}(3);
    p21freqEst(idx) = freqEst{idx}(4); %p21MLerror(idx) = MLerror{idx}(4);
end

f = figure('Position',[500 100 580 640]);
tiledlayout(numExp,4,'TileSpacing','compact','Padding','compact');

panelNames = ["A." "B." "C." "D."];
for idx = 1:numExp
    for idy = 1:4 % plotting the result of the first experiment
    
        nexttile
        [xbar,binCounts] = utilB.bincounts(round(sqrt(length(BayesSeq{idx}))), BayesEst{idx}(idy), BayesError{idx}(idy), BayesPos{idx}(:,idy), BayesPW{idx});
        bar(xbar,binCounts,1,'LineStyle','none','FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',.9) % weighted histogram
        if idy == 1
            xlabel('D_1  (\mum^2/s)')
            ylabel('P(D_1|O)')
            hold on
            xline(D1vec(idx),'--','LineWidth',1.5,'Color',[0,0,0])
            hold on
            xline(freqEst{idx}(1),'-','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980])
            text(-25,102,{panelNames(idx)},'FontSize',12,'unit','pixel'); % text box for indicating panel A,B,C
        elseif idy == 2
            xlabel('D_2  (\mum^2/s)')
            ylabel('P(D_2|O)')
            hold on
            xline(D2vec(idx),'--','LineWidth',1.5,'Color',[0,0,0])
            hold on
            xline(freqEst{idx}(2),'-','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980])
        elseif idy == 3
            xlabel('p_{12}')
            ylabel('P(p_{12}|O)')
            hold on
            xline(p12vec(idx),'--','LineWidth',1.5,'Color',[0,0,0])
            hold on
            xline(freqEst{idx}(3),'-','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980])
        else
            xlabel('p_{21}')
            ylabel('P(p_{21}|O)')
            hold on
            xline(p21vec(idx),'--','LineWidth',1.5,'Color',[0,0,0])
            hold on
            xline(freqEst{idx}(4),'-','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980])
        end
    end
end



