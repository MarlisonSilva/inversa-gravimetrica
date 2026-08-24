function [ gz ] = hypGrav2DPrism( rho, beta, x0, z0, w, h, xobs )
% Local variables
CG = 6.672E-08;
PI = 3.141593;
          
x = xobs - x0;
C = 13.3333*rho*beta^2;        

r1 = sqrt((x + w/2)^2 + z0^2);
r2 = sqrt((x - w/2)^2 + z0^2);        
r3 = sqrt((x - w/2)^2 + (z0 + h)^2);        
r4 = sqrt((x + w/2)^2 + (z0 + h)^2);       

P1 = beta^2 + (x - w/2)^2;
P2 = beta^2 + (x + w/2)^2;

alpha1 = (x - w/2)/P1;
alpha2 = (x + w/2)/P2;

Q1 = (x - w/2)^2 - beta*z0;
Q2 = (x + w/2)^2 - beta*z0;
Q3 = (x - w/2)^2 - beta*(z0 + h);
Q4 = (x + w/2)^2 - beta*(z0 + h);

R = (beta + z0 + h)/(beta + z0);

phi1 = PI/2.0 + atan((x + w/2)/z0);
phi2 = PI/2.0 + atan((x - w/2)/z0);
phi3 = PI/2.0 + atan((x - w/2)/(z0 + h));
phi4 = PI/2.0 + atan((x + w/2)/(z0 + h));

A = R*r2/r3;
B = R*r1/r4;

gz = C*( alpha1*log(A) - alpha2*log(B) + (phi1*Q2/P2 - phi2*Q1/P1)/(beta + z0) + ...
     (phi3*Q3/P1 - phi4*Q4/P2)/(beta + z0 + h) );
end