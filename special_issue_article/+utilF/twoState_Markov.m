
function sVec = twoState_Markov(p12, p21, N)

    % Generates a Markov chain of 1's and 2's.
    % 
    % Input:
    % p12 = transition probability from state 1 to state 2
    % p21 = transition probability from state 2 to state 1
    % N = length of Markov chain
    % 
    % Output:
    % sVec = a row vector [1,1,2,2,2,1,...]
    % 

    pi_1 = p21/(p12+p21); % stationary probability of state 1

    sVec = zeros(1,N); % initialise vector of states
    uVec = rand(1,N); % generate random numbers

    % assign the first state
    if uVec(1) <= pi_1
        sVec(1) = 1; % assign assign the first state (i.e. s1) to be state 1
    else
        sVec(1) = 2; % assign assign the first state (i.e. s1) to be state 2
    end

    % assign remaining states
    desVec1 = uVec <= p12;
    desVec2 = uVec <= p21;
    for i = 2:N
        if sVec(i-1) == 1 % from state 1
            sVec(i) = 2*desVec1(i)+1*~desVec1(i);
        else % from state 2
            sVec(i) = 1*desVec2(i)+2*~desVec2(i);
        end
    end
end
    













