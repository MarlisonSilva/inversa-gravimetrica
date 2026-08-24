clear all;clc;
maxIter=400;            % Maximum number of iterations.
flambda=0.001;		% Marquardt parameter
usrTOL=0.0001;          % Threshold for objective function.  

noiseType=2;
noiseMean=0;
noiseSD=0.1;
noiseSeed=8;

alpha=1.0;		% Total variation regularizing parameter
betaTV=0.001;		% Used to restore differentiability of total variation functional

nprisms=100;            % Number of prisms in the interpretation model.
xminPrisms=0.0;         % Start of discretization along profile.
xmaxPrisms=200.0;       % End of the discretization along profile.
prismDensity=-0.3;	% g/cm3
dcType=1;		% Hyperbolic density contrast 0 - Constant density contrast 1 - Hyperbolic variation with depth.
dcBeta=3.0;             % Decaying factor if dcType = 1.

udist=0;		% 0 - Distances in kilometers, 1 - Distances em meters

P=inv_1D( maxIter,flambda,usrTOL,noiseType,noiseMean,noiseSD,noiseSeed,alpha,betaTV,nprisms,xminPrisms, ...
                         xmaxPrisms,prismDensity,dcType,dcBeta,udist);
P=-1*P;       

clf('reset');
plotarmodelo2D(0,xmaxPrisms,P);