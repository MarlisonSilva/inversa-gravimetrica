&iterations
maxIter=400
flambda=0.001			! Marquardt parameter
usrTOL=0.0001
/
&noiseOptions
noiseType=2
noiseMean=0.0
noiseSD=0.0
noiseSeed=2
/
&regularization
alpha=1.25		! Total variation regularizing parameter
betaTV=0.001	! Used to restore differentiability of total variation functional
/
&model
nprisms=90
xminPrisms=0.0
xmaxPrisms=45
prismDensity=-0.98	! g/cm3
dcType=1			! Hyperbolic density contrast
dcBeta=2.597
/
&units
udist=0			! 0 - Distances in kilometers, 1 - Distances em meters
/