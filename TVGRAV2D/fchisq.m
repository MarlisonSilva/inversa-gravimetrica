function [ Obj, chisqA ] = fchisq(anomObs, YFIT, P, alpha, betaTV, nprismsY, nobsX, nobsY)
%     !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%     ! Objective function
%     !                                 
%     !        ||  o        ||2          
%     ! T(P) = || Y  - F(P) || + alpha*J  (P) 
%     !        ||           ||2         TV
%     !                                 
%     ! Parameters:
%     !   anomObs     Input - Observed data
%     !   YFIT        Input - Data fit for P ( F(P) )
%     !   P           Input - Solution of the inverse problem (prisms thicknesses)
%     !   alpha       Input - Regularization parameter for the total variation functional
%     !   betaTV      Input - Used to restore differentiability of the total variation functional   
%     !   nprismsY    Input - Number of prisms along Y (E-W)
%     !   chisqA      Output - Misfit value normalized by the number of observations.
%     ! Return
%     !    Objective function value
%     !
%     !
%     !//////////////////////////////////////////////////////////////////////////////////////////

chisq = 0;    
k = 0;
w=0;
  for i=1:nobsX
   for j=1:nobsY
    w=j+(i-1)*nobsY;       
    if(anomObs(w,3) ~= -9999.00)
      chisq = chisq + (anomObs(w,3)-YFIT(i,j))^2;					
      k = k + 1;									
    end
   end 
   end 
   
chisqA = sqrt(chisq/(k));
Obj = chisq + alpha*functnTV(P, betaTV, nprismsY);
end