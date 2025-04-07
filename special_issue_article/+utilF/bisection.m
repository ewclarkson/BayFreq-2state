
function [argMinval, nIter] = bisection(f, xL, xR, eps)
%
% Find zeros of a function using the bisection method
%
% Input: 
% f = function handle to function for which we seek f(s) = 0.
% xL = left boundary 
% xR = right boundary
% eps = tolerance
% 
% Output:
% argMinVal = value of argument for which f is zero
% nIter = number of iterations
%
% Code from:
% http://hplgit.github.io/Programming-for-Computations/pub/p4c/._p4c-bootstrap-Matlab028.html
%
    if f(xL)*f(xR) > 0
        fprintf('Error! Function does not have opposite\n');
        fprintf('signs at interval endpoints!')
        return
    end
    xM = (xL + xR)/2.0;
    fM = f(xM);
    iterationCounter = 1;
    while abs(fM) > eps
        left_f = f(xL);
        right_f = f(xR);
        if left_f*fM > 0   % i.e., same sign
            xL = xM;
        else
            xR = xM;
        end
        xM = (xL + xR)/2;
        fM = f(xM);
        iterationCounter = iterationCounter + 1;
    end
    argMinval = xM;
    nIter = iterationCounter;


end

