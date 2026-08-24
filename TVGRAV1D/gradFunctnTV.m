function [ grad ] = gradFunctnTV( P, beta )
%     !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%     ! gradFunctnTV
%     ! Computes the gradient of the total variation functional numerically by 
%     ! simple finite differences approximation.
%     !
%     ! Parameters:
%     !   P       Input - Solution of our inverse problem (prisms thicknesses).
%     !   beta    Input - Parameter used to restore differentiability.
%     !
% 	! Return  
% 	!    Vector, same size as P with the the gradient components.
% 	!
%     !//////////////////////////////////////////////////////////////////////////////

%         !\\\\\\\\\\\\\\\\\\\\\\\\\\
%         ! Start of execution block
%         !//////////////////////////

%tempP=zeros(length(P),1);
grad=zeros(length(P),1);

% for i = 1:length(P)
tempP = P;
% end

 for i = 1:length(P)
   aj = tempP(i);
   dp = 0.01*tempP(i);

  if(dp < 0.001) 
    dp = 0.001; 
  end
            
  tempP(i) = aj+dp;
  a = functnTV(tempP, beta);
  tempP(i) = aj - dp;
  b = functnTV(tempP, beta);
  grad(i) = (a-b)/(2*dp);
  tempP(i) = aj;
end
 
end