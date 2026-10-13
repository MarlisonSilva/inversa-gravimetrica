function [ chisqR, chisqA, tvValue ] = fchisq(GOBS,YFIT,P,alpha,betaTV)
%     !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%     ! fchisq                          
%     ! Objective function
%     !                                 
%     !        ||  o        ||2          
%     ! T(P) = || Y  - F(P) || + alpha*J  (P) 
%     !        ||           ||2         TV
%     !                                 
%     ! Parameters:
%     !   GOBS    Input - Observed data
%     !   YFIT    Input - Data fit for P ( F(P) )
%     !   P       Input - Solution of the inverse problem (prisms thicknesses)
%     !   alpha   Input - Regularization parameter for the total variation functional
%     !   betaTV  Input - Used to restore differentiability of the total variation functional   
%     !   chisqA  Output - Misfit value normalized by the number of observations.
% 	! Return
% 	!    Objective function value
%     !
%     !//////////////////////////////////////////////////////////////////////////////////////////
npts = length(GOBS);    
chisq = 0;
for i = 1:npts
 chisq = chisq + ( GOBS(i)-YFIT(i) )^2;
end

chisqA = sqrt(chisq/npts);
tvValue = functnTV(P, betaTV);
chisqR = chisq + alpha*tvValue;

end