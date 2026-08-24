&iterations
maxIter=400
flambda=0.001			! Marquardt parameter
usrTOL=0.0005
/
&noiseOptions
noiseType=2
noiseMean=0.
noiseSD=0.
noiseSeed=0
/
&regularization
alpha=15.0		! Total variation regularizing parameter
betaTV=0.001		! Used to restore differentiability of total variation functional
/
&model
nprisms=60
xminPrisms=0.0
xmaxPrisms=10.5
prismDensity=-0.63		! g/cm3
dcType=1			! Hyperbolic density contrast
dcBeta=3.0
/
&units
udist=0			! 0 - Distances in kilometers, 1 - Distances em meters
/