function displacements = brownian_displacements_2d(stateVec,sigmaB1,sigmaB2) 
    %
    % Simulation of a 2d, two-state Brownian walk. 
    %
    % Input parameters:
    % stateVec = a state vector with elements = 1 (state 1) or = 2 (state 2)               
    % sigmaB1 = std dev for dx and dy when in state 1
    % sigmaB2 = std dev for dx and dy when in state 2
    %
    % Output:
    % displacements = displacements, sqrt(x(t)^2 + y(t)^2)
    % 

    N = length(stateVec);
    tildeSigmaB = sigmaB1*(2-stateVec) + sigmaB2*(stateVec-1);

    % Generate Brownian displacements and positions 
    Dx = tildeSigmaB.*randn(1,N);   % x-direction
    Dy = tildeSigmaB.*randn(1,N); % y-direction

    % Q = sqrt([Dx(t+dt)]^2 + [Dy(t+dt)]^2)
    displacements = sqrt(Dx.^2 + Dy.^2);

   
end


