function [ grad ] = gradFunctnTV(P, betaTV, nprismsY)
%     ! gradFunctnTV
%     ! Calcula o gradiente do vinculo de variaçao total numericamente
%     ! computes the gradient of the total variation functional numerically by 
%     ! simple finite difference approximation.
%     !
%     ! Parameters:
%     !   P           Input - Solution of our inverse problem (prisms thicknesses).
%     !   betaTV      Input - Parameter used to restore differentiability.
%     !   nprismsY    Input - Number of prisms along Y (E-W)
%     !
%     ! Return  
%     !    Vector, same size as P with the the gradient components.
%     !
%     !//////////////////////////////////////////////////////////////////////////////
tempP=zeros(length(P),1);
grad=zeros(length(P),1);

for k = 1:size(P)
 tempP(k) = P(k);
end

for k = 1:length(tempP)
  aj = tempP(k);
  delta = 0.001*tempP(k);

  if (delta < 0.0001) 
   delta = 0.001;
  end
  tempP(k) = aj+delta;
  td = functnTV(tempP, betaTV, nprismsY);				                 
  yf = td;	   
   
  tempP(k) = aj - delta;		
  td = functnTV(tempP, betaTV, nprismsY);
    
  grad(k) = (yf - td)/(2*delta);

  tempP(k) = aj;
end

end