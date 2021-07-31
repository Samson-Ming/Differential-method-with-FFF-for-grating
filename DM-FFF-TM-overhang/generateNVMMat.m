%% Inputs are toeplitz matrices of sin and cos functions to profile and 
% toeplitz matrices of permittivities, 
%% Outputs are auxiliary matrices A,B,C,D of NVM
function [A,B,C,D] = generateNVMMat(cosM,sinM,etaMn,epsMn)
        I=eye(size(cosM));
        Delta=epsMn-etaMn\I;        
        A=Delta*cosM; %Auxiliary factorization matrices
        B=Delta*sinM;
        D=(epsMn-A)\I;
        C=I;
end




