function [posNew, logLnew, reject] = drawlivepoint2D(pointArr, logLworst, deltaT, dispCell, priorLimits)

    % 
    % Uniformly draw a point obeying a likelihood constraint, by estimating 
    % an ellipsoidal isocontour bounding the live points. Accept a proposed 
    % point if it has a likelihood higher than the constraint. This method 
    % is presented in Shaw et al "Efficient Bayesian inference for 
    % multimodal problems in cosmology" Mon. Not. R. Astron. Soc. 378, 
    % 1365–1370 (2007).
    % 
    % Input:
    % livePoints = array of positions of all live points
    % logLworst = log-likelihood constraint, a lower bound
    % deltaT = sampling time interval
    % dispCell = cell of trajectories of displacements
    % 
    % Output:
    % posNew = position of new point
    % logLnew = likelihood of new point
    % reject = number of rejections of proposed points
    % 
    % Dependencies:
    % EnsLikelihood2D.m
    % LogLikelihood2D.m
    % logsumexp2.m


    [nPoints, nDims] = size(pointArr); % extract #points and #dimensions
    import('utilB.EnsLikelihood2D')
    
    % ------------- approximate a bounding ellipsoid ---------------------
    
    mean_p = mean(pointArr); % estimated mean of points
 
    % Estimate covariance matrix statistically
    app_C = cov(pointArr,1); % compute covariance matrix
    inv_appC = inv(app_C); % inverse of covariance matrix
    
    % Find a bounding ellipsoid
    kVec = zeros(1,nPoints); % helps define an ellipsoid
    for i = 1:nPoints % compute a value in kVec for every point
    
        kVec(i) = (pointArr(i,:)-mean_p)*inv_appC*(pointArr(i,:)-mean_p)';
%         kVec(i) = (pointArr(i,:)-mean_p)*(app_C\((pointArr(i,:)-mean_p)'));
    end
    kScale = max(kVec); % pick the largest k
    kNew = kScale*1.06^2; % slightly enlarge the ellipsoid
    
    % ------------ sample uniformly from bounding ellipsoid --------------
    
    % Compute a transformation matrices from ball to ellipsoid
    [R,D] = eig(app_C); % matrix of eigenvectors and of eigenvalues
    T = sqrt(kNew)*R*sqrt(D)*R'; % transformation matrix from ball to ellipsoid

    % -- test correctness --
    % zArr = zeros(nPoints,nDims);
    % for k = 1:nPoints % sample within D-ball
    %     zt = randn(1,nDims); % random vector of standard gaussians
    %     zt = (rand)^(1/nDims)*zt/norm(zt); % sample within unit ball
    %     zArr(k,:) = zt; % store z
    % end
    % 
    % yVec = zeros(nPoints,nDims);
    % for k = 1:nPoints
    %     yVec(k,:) = (T*zArr(k,:)')'+mean_p; % final random deviates
    % end
    % -- end of test --

    logLnew = logLworst-1;
    reject = -1; % number of rejections
    while logLnew <= logLworst % do until we are inside the region of L>Lworst

        % Generate a random deviate within D-ball
        z = randn(1,nDims); % random vector of standard gaussians
        z = rand^(1/nDims)*z/norm(z); % sample within unit ball
                
        % Apply transformation to D-ball random deviates
        y = (T*z')'+mean_p; % final random deviates
        
        logLnew = EnsLikelihood2D(y, deltaT, dispCell); % likelihood of new point
        % Only accept new point if it is within a uniform prior
        for idy = 1:4
            if priorLimits(idy,1) > y(idy) || y(idy) > priorLimits(idy,2)
                logLnew = logLworst-1;
            end
        end

        reject = reject+1;
    end

    if y(2) > y(1) % D2 > D1

        y = [y(2),y(1),y(4),y(3)]; % enforce D1>D2
    end

    posNew = y; % output
end

% scatter(pointArr(:,1), pointArr(:,2)) 
% hold on
% scatter(yVec(:,1), yVec(:,2))

% scatter(pointArr(:,3), pointArr(:,4))
% hold on
% scatter(yVec(:,3), yVec(:,4))


