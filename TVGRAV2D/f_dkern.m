function [ dkern ] = f_dkern(x, y, z)

R = sqrt(x^2 + y^2 + z^2);	    
phi = atan2(R+x+y,z);
       
 if(x == 0)
  dkern = y*z/(R*(R+x)) + 2*phi + 2*z*( ((1./R) - (1./z^2)*(R+x+y))/(1.+((1./z^2)*(R+x+y)^2)));
 elseif (y == 0)
  dkern = x*z/(R*(R+y)) + 2*phi + 2*z*( ((1./R) - (1./z^2)*(R+x+y))/(1.+((1./z^2)*(R+x+y)^2)));
 elseif (x == 0 &&  y == 0)
  dkern = 2*phi + 2*z*( ((1./R) - (1./z^2)*(R+x+y))/(1.+((1./z^2)*(R+x+y)^2)));
 else
  dkern = x*z/(R*(R+y)) + y*z/(R*(R+x)) + 2*phi + ...
              2*z*( ((1./R) - (1./z^2)*(R+x+y))/(1.+((1./z^2)*(R+x+y)^2)));                        
end       
end