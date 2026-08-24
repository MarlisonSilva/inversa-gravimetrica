&iterations
maxIter=400
flambda=0.0001			! Marquardt parameter
usrTOL=0.0001
/
&noiseOptions
noiseType=2
noiseMean=0.
noiseSD=0.1
noiseSeed=8
/
&regularization
alpha=1.0		! Total variation regularizing parameter
betaTV=0.001		! Used to restore differentiability of total variation functional
/
&model
nprisms=100
xminPrisms=0.0
xmaxPrisms=200.0
prismDensity=-0.3		! g/cm3
dcType=1			! Hyperbolic density contrast
dcBeta=3.0
/
&units
udist=0			! 0 - Distances in kilometers, 1 - Distances em meters
/