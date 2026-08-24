&iterations
maxIter=400
flambda=0.001			! Marquardt parameter
usrTOL=0.0001
/
&noiseOptions
noiseType=2
noiseMean=0
noiseSD=0.0
noiseSeed=5
/
&regularization
alpha=1.25			! Total variation regularizing parameter
betaTV=0.001		! Used to restore differentiability of total variation functional
/
&model
nprisms=60
xminPrisms=0.0
xmaxPrisms=12.0
prismDensity=-0.3		! g/cm3
dcType=0
dcBeta=0
/
&units
udist=0			! 0 - Distances in kilometers, 1 - Distances em meters
/