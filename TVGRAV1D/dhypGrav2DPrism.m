function [ dg ] = dhypGrav2DPrism( rho, beta, x0, z0, w, h, xobs )
CG = 6.672E-08;
PI = 3.141593;

x = xobs - x0;
P1 = beta^2 + (x - w/2)^2;
P2 = beta^2 + (x + w/2)^2;
alpha1 = (x - w/2)/P1;
alpha2 = (x + w/2)/P2;

A = -(-4*x^2+4*x*w-w^2+4*beta*z0+4*beta*h)/((4*z0^2+8*z0*h+4*h^2+4*x^2-4*x*w+w^2)*(beta+z0+h));
B = -(-4*x^2-4*x*w-w^2+4*beta*z0+4*beta*h)/((4*h^2+8*z0*h+4*z0^2+4*x^2+4*x*w+w^2)*(beta+z0+h));

D1 = -(x-w/2)*((x-w/2)^2-beta*(z0+h))/((z0+h)^2*(1.+(x-w/2)^2/(z0+h)^2)*P1*(beta+z0+h));
D2 = -(PI/2+atan((x-w/2)/(z0+h)))*beta/(P1*(beta+z0+h));
D3 = -(PI/2+atan((x-w/2)/(z0+h)))*((x-w/2)^2-beta*(z0+h))/(P1*(beta+z0+h)^2);

D4 = -(x+w/2)*((x+w/2)^2-beta*(z0+h))/((z0+h)^2*(1.+(x+w/2)^2/(z0+h)^2)*P2*(beta+z0+h));
D5 = -(PI/2+atan((x+w/2)/(z0+h)))*beta/(P2*(beta+z0+h));
D6 = -(PI/2+atan((x+w/2)/(z0+h)))*((x+w/2)^2-beta*(z0+h))/(P2*(beta+z0+h)^2);
        
D = (D1 + D2 + D3) - (D4 + D5 + D6);

dg = (13.3333*rho*beta^2)*(alpha1*A - alpha2*B + D);
end