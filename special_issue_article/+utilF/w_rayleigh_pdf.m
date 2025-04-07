function out = w_rayleigh_pdf(xVec, f, th1, th2) 

    % Defines the 2-component rayleigh mixture PDF
    % 
    % Input:
    % xVec = vector of step-sizes, i.e. the support of the PDF
    % f = mixture weight for component 1
    % th1 = scale parameter of component 1
    % th2 = scale parameter of component 2
    % 
    % Output:
    % out = Rayleigh mixture model PDF

    P1 = raylpdf(xVec, th1);
    P2 = raylpdf(xVec, th2);
    
    out = f*P1 + (1-f)*P2;
end
