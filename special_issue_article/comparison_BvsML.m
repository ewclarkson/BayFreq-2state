%% Compare the Bayesian and maximum likelihood methods

% generate and analyse data
% rng(3) % for correspondence between long runs
numExp = 3; % number of data sets to analyse
nTraj = 10;         % number of trajectories
N = 3E2;            % number of displacements
tau = 5E-3;
D1vec = linspace(10,10,numExp);
D2vec = [1,2,3]; %linspace(1,3,numExp);
% sigmaB1vec = sqrt(2*D1vec*tau);
% sigmaB2vec = sqrt(2*D2vec*tau);
p12vec = linspace(0.05,0.05,numExp);          
p21vec = linspace(0.05,0.05,numExp);     
data = cell(1,nTraj);

BayesEst = cell(1,numExp); % hold bayesian parameter estimates
BayesError = cell(1,numExp); % hold error estimates
BayesError2std = cell(1,numExp); % hold 2 standard deviation error estimates
BayesError3std = cell(1,numExp); % hold 2 standard deviation error estimates

NSseq = cell(1,numExp); % hold nested sampling output sequence

MLest = cell(1,numExp); % hold ML parameter estimates
MLerror = cell(1,numExp); % hold ML parameter estimates

for idx = 1:numExp

    D1 = D1vec(idx); D2 = D2vec(idx);
    p12 = p12vec(idx); p21 = p21vec(idx);
    sigmaB1 = sqrt(2*D1*tau); sigmaB2 = sqrt(2*D2*tau);

    % generate a data set
    for idxTrack = 1:nTraj  
        stateVec = utilF.twoState_Markov(p12, p21, N); % random state sequence
        data{idxTrack} = utilF.brownian_displacements_2d(stateVec,sigmaB1,sigmaB2);
    end
    
    % ----------------- do bayesian analysis ------------------------
    D1_max = 20; % D1+0.5*D1; % 20; % % known upper limit of D1
    D1_min = 3; %max(D1-0.5*D1,0.01);% 6; % % known lower limit of D1
    D2_max =  8; %D2+0.5*D2; %5;  % % known upper limit of D2
    D2_min =  0.3; % max(D2-0.5*D2,0.01);%0.3;  % % known lower limit of D2
    p12_max = 0.2; % min(5*p12,1); %  % % known upper limit of p12
    p12_min = 0.005; %min(1/5*p12,0.01);  % % known lower limit of p12
    p21_max = 0.2; %min(5*p21,1); %  % % known upper limit of p21
    p21_min = 0.005; %min(1/5*p21,0.01); % % known lower limit of p21 
    % algorithm parameters
    nLive = 3E2; % number of live points
    StopRatio = 1E-4; % stop criterion for evidence
    % priorLimits = [D1_min,D1_max,D2_min,D2_max,p12_min,p12_max,p21_min,p21_max];
    priorLimits = [D1_min,D1_max;D2_min,D2_max;p12_min,p12_max;p21_min,p21_max];

    tic
    [finalSeq, thetaMLE, logZ] = utilB.nestedsampling2D(nLive, StopRatio, priorLimits, tau, data);
    [logZ_error, thetaBayes, errorBayes] = utilB.bayesianestimate(finalSeq,logZ,nLive,false);
    BayesEst{idx} = thetaBayes; % store parameter estimates
    BayesError{idx} = errorBayes; % store error estimates
    NSseq{idx} = finalSeq;
    toc

    % -------------- do ML-analysis - NOTE: check reasonable values of tweak-parameters beforehand 
    num_MC_steps = 1E4+1;
    burnInTime = 1E3;
    D1_a = 7; D1_b = 14; D2_a = 0.3; D2_b = 3; p12_a = 0.01; p12_b = 0.3; p21_a = 0.01; p22_b = 0.3;
    guessArr = [D1_a, D1_b; D2_a, D2_b; p12_a, p12_b; p21_a, p22_b]; % hold limits on random guess for initial parameter values
    if idx == 1
        scaleArr = [0.9, 0.01, 0.0005, 0.0005]; % set this to achieve acceptance rates of 20-40%
    elseif idx == 2
        scaleArr = [1.1, 0.05, 0.00067, 0.0005];
    elseif idx == 3
        scaleArr = [1.1, 0.17, 0.0007, 0.0006];
    else % if considering more than 3 experiments, i.e. numExp >=3
        scaleArr = [1.1, 0.18, 0.0007, 0.0006];
    end
    
    tic
    % [latest_est, mean_est, D1_est, D2_est, p12_est, p21_est, s1, s2, s3, s4] = MCMC(num_MC_steps, burnInTime, scaleArr, guessArr, tau, data);
    [thetaML, errorUpper, errorLower, s1, s2, s3, s4] = MCMC_v2(num_MC_steps, burnInTime, scaleArr, guessArr, tau, data)
    MLest{idx} = thetaML; % store parameter estimates
    MLerrorUpper{idx} = errorUpper-thetaML; % store upper error estimate
    MLerrorLower{idx} = thetaML-errorUpper; % store error estimates
    % MLerror{idx} = errorLower; % store error estimates
    toc
end



%% Plot 1: keep D1 fixed and vary D2

% extract estimates
D1BayesEst = zeros(1,numExp); D1BayesError = zeros(1,numExp);
D2BayesEst = zeros(1,numExp); D2BayesError = zeros(1,numExp);
p12BayesEst = zeros(1,numExp); p12BayesError = zeros(1,numExp);
p21BayesEst = zeros(1,numExp); p21BayesError = zeros(1,numExp);

D1MLest = zeros(1,numExp); D1MLerror = zeros(1,numExp);
D2MLest = zeros(1,numExp); D2MLerror = zeros(1,numExp);
p12MLest = zeros(1,numExp); p12MLerror = zeros(1,numExp);
p21MLest = zeros(1,numExp); p21MLerror = zeros(1,numExp);

for idx = 1:numExp
    D1BayesEst(idx) = BayesEst{idx}(1); D1BayesError(idx) = BayesError{idx}(1);
    D2BayesEst(idx) = BayesEst{idx}(2); D2BayesError(idx) = BayesError{idx}(2);
    p12BayesEst(idx) = BayesEst{idx}(3); p12BayesError(idx) = BayesError{idx}(3);
    p21BayesEst(idx) = BayesEst{idx}(4); p21BayesError(idx) = BayesError{idx}(4);

    D1MLest(idx) = MLest{idx}(1); D1MLerrorUpper(idx) = MLerrorUpper{idx}(1); D1MLerrorLower(idx) = MLerrorLower{idx}(1);
    D2MLest(idx) = MLest{idx}(2); D2MLerrorUpper(idx) = MLerrorUpper{idx}(2); D2MLerrorLower(idx) = MLerrorLower{idx}(2);
    p12MLest(idx) = MLest{idx}(3); p12MLerrorUpper(idx) = MLerrorUpper{idx}(3); p12MLerrorLower(idx) = MLerrorLower{idx}(3); 
    p21MLest(idx) = MLest{idx}(4); p21MLerrorUpper(idx) = MLerrorUpper{idx}(4); p21MLerrorLower(idx) = MLerrorLower{idx}(4);
end

f = figure('Position',[500 200 580 500]);
% text(0.5,0.5,{"C."},'unit','normalized','FontSize',12); % text box for panel name
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')
nexttile
% figure;
% plot bayesian estimates
plot1 = errorbar(D2BayesEst,D1BayesEst,D1BayesError,D1BayesError,D2BayesError,D2BayesError,...
    'o','Color','#0072BD','MarkerSize',3,'MarkerFaceColor','#0072BD','LineWidth',1.0);
hold on
for idx = 1:numExp % plot estimations of diffusion constants

    finalSeq = NSseq{idx};
    seqLen = length(finalSeq);
    nBins = round(sqrt(seqLen));

    % % include this code in a function instead --------------
    % dataH = zeros(seqLen,4); % all coordinate values
    % w = zeros(1,seqLen); % corresponding weights
    % logLVec = zeros(seqLen,1); % likelihoods
    % 
    % for i = 1:seqLen
    % 
    %     dataH(i,:) = finalSeq(i).pos; % parameter values
    %     w(i) = finalSeq(i).postWt; % posterior weights
    %     logLVec(i) = finalSeq(i).logL; % log-likelihood values
    % end
    % data1 = dataH(:,1); % 1st coordinate
    % data2 = dataH(:,2); % 2nd coordinate
    % data3 = dataH(:,3); % 3rd coordinate
    % data4 = dataH(:,4); % 4th coordinate
    % 
    % [xbar,bc1D] = utilB.bincounts(nBins, D2BayesEst(idx), D2BayesError(idx), data2, w);
    % % bar(xbar,bc1D,1,'LineStyle','none','FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',.9) % weighted histogram
    % % [xbar,binCounts] = utilB.bincounts(nBins, D1BayesEst(idx), D1BayesError(idx), data1, w);
    % % bar(xbar,binCounts,1,'LineStyle','none','FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',.9) % weighted histogram
    % 
    % [~,minIndex] = min(abs(xbar-BayesEst{idx}(2)+BayesError{idx}(2)));
    % lv1 = bc1D(minIndex); 
    % [~,minIndex] = min(abs(xbar-BayesEst{idx}(2)+2*BayesError{idx}(2)));
    % lv2 = bc1D(minIndex); 
    % [~,minIndex] = min(abs(xbar-BayesEst{idx}(2)+3*BayesError{idx}(2)));
    % lv3 = bc1D(minIndex); 
    % % -------------------------------------------------------------------

    nBins2D = round(sqrt(nBins));
    % [xbar1,xbar2,binCounts] = utilB.isocontour(nBins2D, [BayesEst{idx}(1),BayesEst{idx}(2)], [BayesError{idx}(1),BayesError{idx}(2)], data1, data2, w);
    [xbar1,xbar2,binCounts] = utilB.isocontour(finalSeq, nBins2D, [BayesEst{idx}(1),BayesEst{idx}(2)], [BayesError{idx}(1),BayesError{idx}(2)], [1,2]);
    bc = binCounts;
    bcGF = bc;
    % bcGF = imgaussfilt(bc,1.5); % gaussian smoothening
    PmaxGF = max(bcGF,[],'all'); % define maximum
    lv1 = PmaxGF*exp(-1/2);
    lv2 = PmaxGF*exp(-2^2/2); 
    lv3 = PmaxGF*exp(-3^2/2);
    % lv1 = 2.3163;
    % lv2 = 0.5273;
    % lv3 = 0.0332;

    % with filled contours
    contourf(xbar2,xbar1,bcGF, [lv1,lv1],'EdgeColor','none','FaceColor',[0 0.4 1],'FaceAlpha',0.3)
    hold on
    contourf(xbar2,xbar1,bcGF, [lv2,lv2],'EdgeColor','none','FaceColor',[0 0.4 1],'FaceAlpha',0.3)
    hold on
    contourf(xbar2,xbar1,bcGF, [lv3,lv3],'EdgeColor','none','FaceColor',[0 0.4 1],'FaceAlpha',0.3)
end

plot0 = plot(D2vec,D1vec,'o','MarkerSize',4,'MarkerFaceColor','#D95319','MarkerEdgeColor','none');
hold on

% plot ML estimates
plot2 = errorbar(D2MLest,D1MLest,D1MLerrorLower,D1MLerrorUpper,D2MLerrorLower,D2MLerrorUpper,...
    'o','Color','#77AC30','MarkerSize',3,'MarkerFaceColor','#77AC30','LineWidth',1.0);

xlabel('D_2  (\mum^2/s)')
ylabel('D_1  (\mum^2/s)')
% legend('Ground truth','Bayesian estimate','ML estimate')
legend([plot0,plot1,plot2],{'Ground truth','ns-HMM','mh-HMM'},'AutoUpdate','off','Location','northwest')
grid on
text(-0.11,0.97,{"A."},'unit','normalized','FontSize',12); % text box for panel name

% Plot p12 vs D2
nexttile
% plot bayesian estimates
plot1 = errorbar(D2BayesEst,p12BayesEst,p12BayesError,p12BayesError,D2BayesError,D2BayesError,...
'o','Color','#0072BD','MarkerSize',3,'MarkerFaceColor','#0072BD','LineWidth',1.0);

hold on
for idx = 1:numExp % plot estimations of diffusion constants

    finalSeq = NSseq{idx};
    seqLen = length(finalSeq);
    nBins = round(sqrt(seqLen));

    % % include this code in a function instead --------------
    % dataH = zeros(seqLen,4); % all coordinate values
    % w = zeros(1,seqLen); % corresponding weights
    % logLVec = zeros(seqLen,1); % likelihoods
    % 
    % for i = 1:seqLen
    % 
    %     dataH(i,:) = finalSeq(i).pos; % parameter values
    %     w(i) = finalSeq(i).postWt; % posterior weights
    %     logLVec(i) = finalSeq(i).logL; % log-likelihood values
    % end
    % data1 = dataH(:,1); % 1st coordinate
    % data2 = dataH(:,2); % 2nd coordinate
    % data3 = dataH(:,3); % 3rd coordinate
    % data4 = dataH(:,4); % 4th coordinate
    % 
    % [xbar,bc1D] = utilB.bincounts(nBins, D2BayesEst(idx), D2BayesError(idx), data2, w);
    % % bar(xbar,bc1D,1,'LineStyle','none','FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',.9) % weighted histogram
    % % [xbar,binCounts] = utilB.bincounts(nBins, D1BayesEst(idx), D1BayesError(idx), data1, w);
    % % bar(xbar,binCounts,1,'LineStyle','none','FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',.9) % weighted histogram
    % 
    % [~,minIndex] = min(abs(xbar-BayesEst{idx}(2)-BayesError{idx}(2)));
    % lv1 = bc1D(minIndex);
    % % lv1 = 100;
    % [~,minIndex] = min(abs(xbar-BayesEst{idx}(2)-2*BayesError{idx}(2)));
    % lv2 = bc1D(minIndex); 
    % [~,minIndex] = min(abs(xbar-BayesEst{idx}(2)-3*BayesError{idx}(2)));
    % lv3 = bc1D(minIndex); 
    % % ------------------------------------------------------------------

    nBins2D = round(sqrt(nBins));
    [xbar1,xbar2,binCounts] = utilB.isocontour(finalSeq, nBins2D, [BayesEst{idx}(3),BayesEst{idx}(2)], [BayesError{idx}(3),BayesError{idx}(2)], [3,2]);
    bc = binCounts;
    bcGF = bc;
    % bcGF = imgaussfilt(bc,1.5); % gaussian smoothening
    PmaxGF = max(bcGF,[],'all'); % define maximum
    lv1 = PmaxGF*exp(-1/2);
    lv2 = PmaxGF*exp(-2^2/2); 
    lv3 = PmaxGF*exp(-3^2/2);
    % lv1 = 0.005;
    % lv2 = 0.0012;
    % lv3 = 0.0001;

    % [X,Y] = meshgrid(xbar1,xbar2);
    % bar3(xbar1,binCounts,1)

    % with filled contours
    contourf(xbar2,xbar1,bcGF, [lv1,lv1],'EdgeColor','none','FaceColor',[0 0.4 1],'FaceAlpha',0.3)
    hold on
    contourf(xbar2,xbar1,bcGF, [lv2,lv2],'EdgeColor','none','FaceColor',[0 0.4 1],'FaceAlpha',0.3)
    hold on
    contourf(xbar2,xbar1,bcGF, [lv3,lv3],'EdgeColor','none','FaceColor',[0 0.4 1],'FaceAlpha',0.3)
end

plot0 = plot(D2vec,p12vec,'o','MarkerSize',4,'MarkerFaceColor','#D95319','MarkerEdgeColor','none');
hold on

% plot ML estimates
plot2 = errorbar(D2MLest,p12MLest,p12MLerrorLower,p12MLerrorUpper,D2MLerrorLower,D2MLerrorUpper,...
    'o','Color','#77AC30','MarkerSize',3,'MarkerFaceColor','#77AC30','LineWidth',1.0);

xlabel('D_2  (\mum^2/s)')
ylabel('p_{12}')
grid on
text(-0.11,0.97,{"B."},'unit','normalized','FontSize',12); % text box for panel name



%%
% testing 1D histograms

% dataH = zeros(seqLen,4); % all coordinate values
% w = zeros(1,seqLen); % corresponding weights
% logLVec = zeros(seqLen,1); % likelihoods
% 
% for i = 1:seqLen
% 
%     dataH(i,:) = finalSeq(i).pos; % parameter values
%     w(i) = finalSeq(i).postWt; % posterior weights
%     logLVec(i) = finalSeq(i).logL; % log-likelihood values
% end
% data1 = dataH(:,1); % 1st coordinate
% data2 = dataH(:,2); % 2nd coordinate
% data3 = dataH(:,3); % 3rd coordinate
% data4 = dataH(:,4); % 4th coordinate
% 
% [xbar,bc1D] = utilB.bincounts(nBins, D2BayesEst(idx), D2BayesError(idx), data2, w);
% bar(xbar,bc1D,1,'LineStyle','none','FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',.9) % weighted histogram
% % [xbar,binCounts] = utilB.bincounts(nBins, D1BayesEst(idx), D1BayesError(idx), data1, w);
% % bar(xbar,binCounts,1,'LineStyle','none','FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',.9) % weighted histogram
% 
% [~,minIndex] = min(abs(xbar-BayesEst{idx}(2)+BayesError{idx}(2)));
% lv1 = bc1D(minIndex); 
% [~,minIndex] = min(abs(xbar-BayesEst{idx}(2)+2*BayesError{idx}(2)));
% lv2 = bc1D(minIndex); 
% [~,minIndex] = min(abs(xbar-BayesEst{idx}(2)+3*BayesError{idx}(2)));
% lv3 = bc1D(minIndex); 



%%
% tc = randn(1,400);
% plot(tc);
% [f,xi] = ksdensity(tc); 
% plot(xi,f);


%% plot p21 vs D2
% 
% figure;
% plot1 = errorbar(D2BayesEst,p21BayesEst,p21BayesError,p21BayesError,D2BayesError,D2BayesError,...
% 'o','Color','#0072BD','MarkerSize',3,'MarkerFaceColor','#0072BD','LineWidth',1.0);
% 
% hold on
% for idx = 1:numExp % plot estimations of diffusion constants
% 
%     finalSeq = NSseq{idx};
%     seqLen = length(finalSeq);
%     nBins = round(sqrt(seqLen));
% 
%     nBins2D = round(sqrt(nBins));
%     [xbar1,xbar2,binCounts] = utilB.isocontour(finalSeq, nBins2D, [BayesEst{idx}(4),BayesEst{idx}(2)], [BayesError{idx}(4),BayesError{idx}(2)], [4,2]);
%     bc = binCounts;
%     bcGF = bc;
%     % bcGF = imgaussfilt(bc,1.0); % gaussian smoothening
%     PmaxGF = max(bcGF,[],'all'); % define maximum
%     lv1 = PmaxGF*exp(-1/2);
%     lv2 = PmaxGF*exp(-2/2); 
%     lv3 = PmaxGF*exp(-3/2);
%     % lv1 = 0.005;
%     % lv2 = 0.0012;
%     % lv3 = 0.0001;
% 
%     % with filled contours
%     contourf(xbar2,xbar1,bcGF, [lv1,lv1],'EdgeColor','none','FaceColor',[0 0.3 1],'FaceAlpha',0.3)
%     hold on
%     contourf(xbar2,xbar1,bcGF, [lv2,lv2],'EdgeColor','none','FaceColor',[0 0.3 1],'FaceAlpha',0.3)
%     hold on
%     contourf(xbar2,xbar1,bcGF, [lv3,lv3],'EdgeColor','none','FaceColor',[0 0.3 1],'FaceAlpha',0.3)
% end
% 
% plot0 = plot(D2vec,p21vec,'o','MarkerSize',4,'MarkerFaceColor','#D95319','MarkerEdgeColor','none');
% hold on

% % plot ML estimates
% plot2 = errorbar(D2MLest,p21MLest,p21MLerror,p21MLerror,D2MLerror,D2MLerror,...
%     'o','Color','#77AC30','MarkerSize',3,'MarkerFaceColor','#77AC30','LineWidth',1.0);
% 
% xlabel('D_2  (\mum^2/s)')
% ylabel('p_{21}')
% grid on
% text(-0.11,0.97,{"B."},'unit','normalized','FontSize',12); % text box for panel name
% 


%%


    %     tEl = 0:0.001:2*pi;
    % --------------- Diffusion constants --------------------
%     % first confidence bound, about 68%
%     xEl = D2BayesEst(idx) + D2BayesError(idx)*sin(tEl);
%     yEl = D1BayesEst(idx) + D1BayesError(idx)*cos(tEl);
%     fill(xEl,yEl,[0,0.3,1],'FaceAlpha',0.2,'EdgeColor','none')
%     hold on
%     % second confidence bound, about 95%
%     xEl = D2BayesEst(idx) + 2*D2BayesError(idx)*sin(tEl);
%     yEl = D1BayesEst(idx) + 2*D1BayesError(idx)*cos(tEl);
%     fill(xEl,yEl,[0,0.3,1],'FaceAlpha',0.3,'EdgeColor','none')
%     hold on
%     % third confidence bound, about 100%
%     xEl = D2BayesEst(idx) + 3*D2BayesError(idx)*sin(tEl);
%     yEl = D1BayesEst(idx) + 3*D1BayesError(idx)*cos(tEl);
%     fill(xEl,yEl,[0,0.3,1],'FaceAlpha',0.3,'EdgeColor','none')
%     hold on


%%
    % % extract data1 and data2, i.e. all values of D1 and D2
    % dataT = zeros(seqLen,4); % all coordinate values
    % w = zeros(1,seqLen); % corresponding weights
    % logLVec = zeros(seqLen,1); % likelihoods
    % for i = 1:seqLen
    % 
    %     dataT(i,:) = finalSeq(i).pos; % parameter values
    %     w(i) = finalSeq(i).postWt; % posterior weights
    %     logLVec(i) = finalSeq(i).logL; % log-likelihood values
    % end
    % data1 = dataT(:,1); % 1st coordinate
    % data2 = dataT(:,2); % 2nd coordinate

        
        % plot Bayesian results for day2,...,day10
% nDays = 5;
% for idd = 1:nDays
%     for idx = 1:numExp
%     
%         D1 = D1vec(idx); D2 = D2vec(idx);
%         p12 = p12vec(idx); p21 = p21vec(idx);
%         sigmaB1 = sqrt(2*D1*tau); sigmaB2 = sqrt(2*D2*tau);
%     
%         % generate a data set
%         dataMoreDays = cell(1,nTraj);
%         for idxTrack = 1:nTraj  
%             stateVec = utilF.twoState_Markov(p12, p21, N); % random state sequence
%             dataMoreDays{idxTrack} = utilF.brownian_displacements_2d(stateVec,sigmaB1,sigmaB2);
%         end
%         
%         % ----------------- do bayesian analysis ------------------------
%         D1_max = D1+0.5*D1; % known upper limit of D1
%         D1_min = max(D1-0.5*D1,0.01); % known lower limit of D1
%         D2_max = D2+0.5*D2; % known upper limit of D2
%         D2_min = max(D2-0.5*D2,0.01); % known lower limit of D2
%         p12_max = min(5*p12,1); % known upper limit of p12
%         p12_min = min(1/5*p12,0.01); % known lower limit of p12
%         p21_max = min(5*p21,1); % known upper limit of p21
%         p21_min = min(1/5*p21,0.01); % known lower limit of p21 
%         % algorithm parameters
%         nLive = 100; % number of live points
%         StopRatio = 1E-4; % stop criterion for evidence
%         priorLimits = [D1_min,D1_max,D2_min,D2_max,p12_min,p12_max,p21_min,p21_max];
%     
%         [finalSeq, thetaMLE, logZ] = utilB.nestedsampling2D(nLive, StopRatio, priorLimits, tau, dataMoreDays);
%         [logZ_error, thetaBayes, errorBayes] = utilB.bayesianestimate(finalSeq,logZ,nLive,false);
%     %     BayesEst{idx} = thetaBayes; % store parameter estimates
%     %     BayesError{idx} = errorBayes; % store error estimates
%     %     
%         hold on
%         plot(thetaBayes(2),thetaBayes(1),'x','MarkerFaceColor','black', 'MarkerEdgeColor','black')
%     
%         % -------------- do ML-analysis - NOTE: check reasonable values of tweak-parameters beforehand 
%     %     num_MC_steps = 1E4+1;
%     %     burnInTime = 1E3;
%     %     scaleArr = [1.4, 0.02, 0.002, 0.002]; % set this to achieve acceptance rates of 20-40%
%     %     D1_a = 7; D1_b = 14; D2_a = 0.3; D2_b = 3; p12_a = 0.01; p12_b = 0.3; p21_a = 0.01; p22_b = 0.3;
%     %     guessArr = [D1_a, D1_b; D2_a, D2_b; p12_a, p12_b; p21_a, p22_b]; % hold limits on random guess for initial parameter values  
%         
%     
%     %     [thetaML, errorML, s1, s2, s3, s4] = MCMC_v2(num_MC_steps, burnInTime, scaleArr, guessArr, tau, data)
%     %     MLest{idx} = thetaML; % store parameter estimates
%     %     MLerror{idx} = errorML; % store error estimates
%     
%     end
% end


% nexttile
% % plot transition probabilities
% for idx = 1:numExp
%     % ------------- Transition probabilities ------------------
%     xEl = p21BayesEst(idx) + p21BayesError(idx)*sin(tEl);
%     yEl = p12BayesEst(idx) + p12BayesError(idx)*cos(tEl);
%     fill(xEl,yEl,[0,0.3,1],'FaceAlpha',0.2,'EdgeColor','none')
%     hold on
%     % second confidence bound, about 95%
%     xEl = p21BayesEst(idx) + 2*p21BayesError(idx)*sin(tEl);
%     yEl = p12BayesEst(idx) + 2*p12BayesError(idx)*cos(tEl);
%     fill(xEl,yEl,[0,0.3,1],'FaceAlpha',0.3,'EdgeColor','none')
%     hold on
%     % third confidence bound, about 100%
%     xEl = p21BayesEst(idx) + 3*p21BayesError(idx)*sin(tEl);
%     yEl = p12BayesEst(idx) + 3*p12BayesError(idx)*cos(tEl);
%     fill(xEl,yEl,[0,0.3,1],'FaceAlpha',0.3,'EdgeColor','none')
%     hold on
% end


% xlim([0.1,5])
% ylim([0,2*D1])



%% Plot 2: keep p12 fixed and vary p21


% construct a matrix/mesh with correct bin centra

% figure;
% [xbar,binCounts] = utilB.bincounts(nBins, D1BayesEst(1), D1BayesError(1), data1, w);
% bar(xbar,binCounts,1)
% 
% figure;
% [xbar,binCounts] = utilB.bincounts(nBins, D2BayesEst(1), D2BayesError(1), data2, w);
% bar(xbar,binCounts,1)










