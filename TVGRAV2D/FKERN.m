function [ Kern ] = FKERN(X,Y,Z)
%     !C FUNCTION FKERN COMPUTE THE WORK FUNCTION
%     !C X*LOG(R+Y)+Y*LOG(R+X)+2*Z*ATAN((R+X+Y)/Z)
R=sqrt(X*X+Y*Y+Z*Z);
CREAL=Z;
CIMAG=R+X+Y;
PHI=atan2(CIMAG,CREAL);
    
    if (X == 0) 
       Kern= Y*log(R+X)+2*Z*PHI;
    elseif (Y == 0)
       Kern= X*log(R+Y)+2*Z*PHI;
    elseif (X == 0 && Y ==0)
       Kern= 2*Z*PHI;
    else
       Kern=X*log(R+Y)+Y*log(R+X)+2*Z*PHI;
    end
    
return;

end