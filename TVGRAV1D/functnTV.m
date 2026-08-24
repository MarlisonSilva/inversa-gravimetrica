function [ TV ] = functnTV( P, beta )
%     !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%     ! functnTV
%     ! Total Variation functional, given by the equation below.
% 	!
%     !  m - 1
%     !  ----    ______________________
%     !   \ ´   |            2       2
%     !   /   - | |P    - P |  + beta
%     !  /__\  \|   i+1    i
%     !  i = 1
%     !
% 	! Parameters
% 	!	P	    Input - Solution of our inverse problem (prisms thicknesses).
%     !   beta    Input - Parameter used to restore differentiability.
%     !
% 	! Return
% 	!   Total variation of the solution P.
%     !
%     !////////////////////////////////////////////////////////////////////////////

TV = 0.0;
for i = 1:size(P) - 1
 dp = P(i+1) - P(i);
 TV = TV + sqrt( dp^2 + beta^2 );
end 

end