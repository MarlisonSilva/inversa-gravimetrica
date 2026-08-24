function Fmod1= modGrav1dTOTAL( P, XOBS, xminPrisms, prismWidth, prismDensity )  

for iobs=1:length(XOBS)
    Fmod1(iobs)= modGrav1d2( P, XOBS, iobs, xminPrisms, prismWidth, prismDensity ) ;
end

%ANOM=Fmod1;
end