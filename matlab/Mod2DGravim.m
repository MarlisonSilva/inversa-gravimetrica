function Dgrav=Mod2DGravim(m)
nprisms=60;	%! Number of prisms in the interpretation model.
xminPrisms=0;	%! Start of discretization along profile.
xmaxPrisms=60;	%! End of the discretization along profile.
prismWidth=1;%! Prisms horizontal width, the same for all prisms.
prismDensity=-0.25;	%! g/cm3! Density contrast of the prisms, the same for all prisms.
nobs = 60; % Mudar se necessário
xobs = 1:60;


Dgrav = modGrav1dTOTAL( m, xobs, xminPrisms, prismWidth, prismDensity ) ;



end