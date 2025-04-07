%% Compare the Bayesian and frequentist approaches

%
numExp = 4; % number of data sets to analyse
nTraj = 100;         % number of trajectories
N = 3E2;            % number of displacements
tau = 5E-3;
D1vec = [10,10,10,10]; %linspace(10,10,numExp);
D2vec = [1,3,1,3]; %linspace(1,4,numExp);
% sigmaB1vec = sqrt(2*D1vec*tau);
% sigmaB2vec = sqrt(2*D2vec*tau);
p12vec = [0.01 0.01 0.05 0.05]; %linspace(0.01,0.01,numExp);          
p21vec = [0.01 0.01 0.05 0.05]; %linspace(0.01,0.01,numExp);
% fTrueVec = p21vec./(p12vec+p21vec);
data = cell(1,nTraj);

BayesEst = cell(1,numExp); % hold bayesian parameter estimates
BayesError = cell(1,numExp); % hold bayesian error estimates
BayesSeq = cell(1,numExp); % full NS output
freqEst = cell(1,numExp); % hold frequentist parameter estimates
% freqError = cell(1,numExp); % hold frequentist error estimates % NOTE: no error estimations available

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
    % [estInitial,~] = utilF.Rayleigh_2state(D1, D2, fTrue, tau, data); % NOTE: how to pick initial guess?
    [estInitial,~] = utilF.Rayleigh_2state((D1_min+D1_max)/2, (D2_min+D2_max)/2, (p21_max/(p21_max+p12_min)+p21_min/(p21_min+p12_max))/2, tau, data); 
    freq_fEst = estInitial(1); freq_D1Est = estInitial(2); freq_D2Est = estInitial(3); % estimates of f,D1,D2
    % [freq_p12Est, freq_p21Est] = utilF.method_of_moments_multi(data, freq_fEst, sqrt(2*freq_D1Est*tau), sqrt(2*freq_D2Est*tau));
   [freq_p12Est, freq_p21Est] = utilF.wrapperFreq(data, freq_fEst, sqrt(2*freq_D1Est*tau), sqrt(2*freq_D2Est*tau));

    freqEst{idx} = [freq_D1Est,freq_D2Est,freq_p12Est,freq_p21Est]; % store parameter estimates
%     freqError{idx} = ; % store error estimates
    toc
end


%%

% extract estimates
% D1BayesEst = zeros(1,numExp); D1BayesError = zeros(1,numExp);
% D2BayesEst = zeros(1,numExp); D2BayesError = zeros(1,numExp);
% p12BayesEst = zeros(1,numExp); p12BayesError = zeros(1,numExp);
% p21BayesEst = zeros(1,numExp); p21BayesError = zeros(1,numExp);

BayesPos = cell(1,numExp); % holds all positions
BayesPW = cell(1,numExp); % holds all posterior weights

% BayesPos = [];
% BayesPW = [];

D1freqEst = zeros(1,numExp); %D1MLerror = zeros(1,numExp);
D2freqEst = zeros(1,numExp); %D2MLerror = zeros(1,numExp);
p12freqEst = zeros(1,numExp); %p12MLerror = zeros(1,numExp);
p21freqEst = zeros(1,numExp); %p21MLerror = zeros(1,numExp);

for idx = 1:numExp
%     D1BayesEst(idx) = BayesEst{idx}(1); D1BayesError(idx) = BayesError{idx}(1);
%     D2BayesEst(idx) = BayesEst{idx}(2); D2BayesError(idx) = BayesError{idx}(2);
%     p12BayesEst(idx) = BayesEst{idx}(3); p12BayesError(idx) = BayesError{idx}(3);
%     p21BayesEst(idx) = BayesEst{idx}(4); p21BayesError(idx) = BayesError{idx}(4);
%     BayesPos{idx} = 
%     BayesPW{idx} =
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

% f = figure('Position',[1000 500 600 330]);
f = figure('Position',[500 100 580 640]);
% f = figure;
tiledlayout(numExp,4,'TileSpacing','compact','Padding','compact');
% tiledlayout(2,2)
% str = {{'A simple plot','from 1 to 10'},'y = x'};
% text([0 0.8],[7 7],str,'unit','pixels')

% idx = 1;
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
            % text(-0.4,1.1,{panelNames(idx)},'unit','normalized','FontSize',12); % text box for indicating panel A,B,C
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

% str = {'A.'};
% t = text(-950,1050,str,'unit','pixels');
% t = text(-4.3,2.2,str,'unit','normalized','FontSize',12);
% t = text(-900,1200,str,'unit','pixels','FontSize',12);
% t = text(0,0,str,'unit','normalized','FontSize',12,'Parent',f);


%%

% [p12DoublyCorrected, p21DoublyCorrected] = utilF.wrapperFreq(data, freq_fEst, sqrt(2*freq_D1Est*tau), sqrt(2*freq_D2Est*tau))

%% Plot step-size distribution with weighted rayleigh and components

% 
%     % ----------------- extract step sizes ----------------
% 
%     StepSizes = []; % initialise holder for all estimated diffusion constants
% 
%     for k = 1:nTraj % do for every track
%         for j = 1:length(data{k}) % do for every displacement
% 
%             dist = data{k}(j); % displacement
%             StepSizes(end+1) = dist; % store displacement
%         end
%     end
% 
%     % --------------- compute a maximum-likelihood estimate -----------------
% 
%     in_guess = [fTrue, sqrt(2*tau*D1), sqrt(2*tau*D2)]; % f, th1, th2
%     mypdf = @(data,f,th1,th2) utilF.w_rayleigh_pdf(data, f, th1, th2); % define function
% 
%     est = mle(StepSizes, 'pdf', mypdf, 'start', in_guess);  
% %     acov = mlecov(est,StepSizes,'pdf',mypdf);
% %     errorEst = sqrt(diag(acov));  % estimated errors for each of the parameters
%     errorEst = NaN;
% 
%     % plotting
%     xVec = linspace(0, max(StepSizes), 1000); % create x-vector for plotting
%     histogram(StepSizes, 'Normalization', 'pdf', 'FaceColor',[1 0 0], 'EdgeColor', [1 0 0], 'LineStyle','none');
% %     xlabel('Step size (\mum^2)')
% %     ylabel('Normalised frequency')
%     hold on
%     plot(xVec, utilF.w_rayleigh_pdf(xVec, fTrue, sqrt(2*D1*tau), sqrt(2*D2*tau)), 'linewidth', 2.0, 'Color','black')
%     hold on
%     plot(xVec, fTrue*raylpdf(xVec,sqrt(2*D1*tau)),'Linewidth', 2.0, 'Color','black','LineStyle','--')
%     hold on
%     plot(xVec, (1-fTrue)*raylpdf(xVec,sqrt(2*D2*tau)),'Linewidth', 2.0, 'Color','black','LineStyle','--')
%     xticklabels([])
%     yticklabels([])
% %     box off


%% Plot 1: keep D1 fixed and vary D2 (p12,p21 fixed and not visible)
% 
% % extract estimates
% D1BayesEst = zeros(1,numExp); D1BayesError = zeros(1,numExp);
% D2BayesEst = zeros(1,numExp); D2BayesError = zeros(1,numExp);
% p12BayesEst = zeros(1,numExp); p12BayesError = zeros(1,numExp);
% p21BayesEst = zeros(1,numExp); p21BayesError = zeros(1,numExp);
% 
% D1freqEst = zeros(1,numExp); %D1MLerror = zeros(1,numExp);
% D2freqEst = zeros(1,numExp); %D2MLerror = zeros(1,numExp);
% p12freqEst = zeros(1,numExp); %p12MLerror = zeros(1,numExp);
% p21freqEst = zeros(1,numExp); %p21MLerror = zeros(1,numExp);
% 
% for idx = 1:numExp
%     D1BayesEst(idx) = BayesEst{idx}(1); D1BayesError(idx) = BayesError{idx}(1);
%     D2BayesEst(idx) = BayesEst{idx}(2); D2BayesError(idx) = BayesError{idx}(2);
%     p12BayesEst(idx) = BayesEst{idx}(3); p12BayesError(idx) = BayesError{idx}(3);
%     p21BayesEst(idx) = BayesEst{idx}(4); p21BayesError(idx) = BayesError{idx}(4);
% 
%     D1freqEst(idx) = freqEst{idx}(1); %D1MLerror(idx) = MLerror{idx}(1);
%     D2freqEst(idx) = freqEst{idx}(2); %D2MLerror(idx) = MLerror{idx}(2);
%     p12freqEst(idx) = freqEst{idx}(3); %p12MLerror(idx) = MLerror{idx}(3);
%     p21freqEst(idx) = freqEst{idx}(4); %p21MLerror(idx) = MLerror{idx}(4);
% end
% 
% % plot ground truth quantities
% % D2GT = linspace(D2,4*D2,4);
% % D1GT = D1*ones(1,4);
% plot(D2vec,D1vec,'o','MarkerFaceColor','#D95319','MarkerEdgeColor','none')
% hold on
% 
% % plot bayesian estimates
% errorbar(D2BayesEst,D1BayesEst,D1BayesError,D1BayesError,D2BayesError,D2BayesError,...
%          'o','MarkerFaceColor','#0072BD','MarkerEdgeColor','none','Color','#0072BD')
% hold on
% 
% % plot frequentist estimates
% errorbar(D2freqEst,D1freqEst,[],[],[],[],'o','MarkerFaceColor','#77AC30','MarkerEdgeColor','none')
% 
% grid on
% % xlim([0.1,5])
% % ylim([0,2*D1])
% xlabel('D_2')
% ylabel('D_1')
% legend('Ground truth','Bayesian estimate','frequentist estimate')
% 



