&iterations
maxIter=400
flambda=0.001			! Marquardt parameter
usrTOL=0.0005
/
&noiseOptions
noiseType=2
noiseMean=0.
noiseSD=0.1
noiseSeed=2
/
&regularization
alpha=1.25		! Total variation regularizing parameter
betaTV=0.001	! Used to restore differentiability of total variation functional
/
&model
nprisms=120
xminPrisms=0.0
xmaxPrisms=60.0
prismDensity=-0.24	! g/cm3
dcType=0			
dcBeta=0
/
&units
udist=0			! 0 - Distances in kilometers, 1 - Distances em meters
/