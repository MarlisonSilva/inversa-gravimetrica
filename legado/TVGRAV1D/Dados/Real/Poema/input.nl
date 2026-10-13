&iterations
maxIter=200
flambda=0.001		! Marquardt parameter
usrTOL=0.0001
/
&noiseOptions
noiseType=2
noiseMean=0
noiseSD=0.0
noiseSeed=1
/
&regularization
alpha=0.001		! Total variation regularizing parameter
betaTV=0.0001		! Used to restore differentiability of total variation functional
/
&model
nprisms=45
xminPrisms=0.0
xmaxPrisms=22.5
prismDensity=-2.3	! g/cm3
dcType=0
dcBeta=0
/
&units
udist=1			! 0 - Distances in kilometers, 1 - Distances em meters
/