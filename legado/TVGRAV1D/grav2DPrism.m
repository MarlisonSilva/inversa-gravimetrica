function [ gz ] = grav2DPrism( den, x0, z0,w, h, xobs )
% !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% 	! grav2DPrism
% 	! Computes the gravity field for a 2D prism with constant density contrast.
% 	! Ref.: Telford, Geldart, Sheriff. Applied Geophysics p. 73, eq. 2.49c.
%     !
%     ! Parameters:
% 	!  den  Input - Density contrast (g/cc)
% 	!  x0   Input - Horizontal centre of the prism
% 	!  z0   Input - Top of the prism
% 	!  w    Input - Prism width
% 	!  h    Input - Prism thickness
% 	!  xobs Input - Observational point where to compute the field
% 	!
% 	! Return:
%     !   Value of the field:
% 	!	A) MILIGAL If density contrast in g/cc and distances in km.
% 	!	B) MICROGAL If density contrast in g/cc and distances in m.
% 	!
% 	! PS. Its assumed z = 0 for observations
% !///////////////////////////////////////////////////////////////////////////////
CG = 6.672e-08;

cd = den*1000;
x = xobs - x0 + w/2.0;

r1 = sqrt(x^2 + z0^2);
r2 = sqrt(x^2 + (h+z0)^2);
r3 = sqrt((x-w)^2 + z0^2);
r4 = sqrt((x-w)^2 + (h+z0)^2);
        
theta1 = atan(x/z0);
theta2 = atan(x/(z0+h));
theta3 = atan((x-w)/z0);
theta4 = atan((x-w)/(z0+h));

gz = 2*10^5*CG*cd*(x*log((r2*r3)/(r1*r4)) + w*log(r4/r3) - h*(theta4-theta2) + z0*(theta3-theta1));
end