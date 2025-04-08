function threshDispl = determine_optimal_thresh_displ(sigmaB1,sigmaB2) 
    %
    % Determine the an "optimal" threshold for the displacements
    % for a two-state Brownian motion in 2d. The threshold is determined
    % so that the 1-CDF_state1(s) = CDF_state2(s) 
    %
    % Input:       
    % sigmaB1 = std dev for dx and dy when in state 1
    % sigmaB2 = std dev for dx and dy when in state 2
    %
    % Output:
    % threshDispl = "optimal" threshold for the displacements.
    %
    % Comment:  
    % The displacements = s = sqrt(x(t)^2 + y(t)^2) 
    % are described by Rayleigh PDF:
    %     PDF(s) = (s/w^2)*exp(-s^2/(2*w^2))
    % where 
    % w = sigmaB1 (when in state 1) and w=sigmaB2 (when in state 2)
    % The cumulative distribution function is:
    %       CDF(s) = 1-exp(-s^2/(2*w^2))
    %                 
    % Dependencies:func_two_rayleighs
    %

    sigmaMin = min(sigmaB1,sigmaB2);
    sigmaMax = max(sigmaB1,sigmaB2);   
    xL = 0;
    xR = 3*sigmaMax;  
                                
    f = @(s) utilF.func_two_rayleighs(s, sigmaMin, sigmaMax);   
    threshDispl = fzero(f,[xL,xR]); % root-finding

    %
    %   % Plot CDFs
    %     s=0:0.01:10;  % displacements   
    %     figure
    %     y1 = 1 - raylcdf(s,sigmaMin);
    %     plot(s,y1,'b-','linewidth',3);   hold on;
    %     y2 = raylcdf(s,sigmaMax);
    %     plot(s,y2,'r-','linewidth',3)
    %     y =[0 1];
    %     x = threshDispl*[1 1];
    %     plot(x,y,'--','linewidth',2);
    %     
    
end


