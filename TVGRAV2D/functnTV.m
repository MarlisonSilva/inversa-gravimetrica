function [ fTV ] = functnTV(P, betaTV, nprismsY)
% !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% ! functnTV
% ! Total Variation functional in 2D
% !
% ! Parameters
% !	P	        Input - Solution of our inverse problem (prisms thicknesses).
% !   betaTV      Input - Parameter used to restore differentiability.
% !   nprismsY    Input - Number of prisms along Y (E-W)
% !
% ! Return
% !   Total Variation of the solution P.
% !
% !
% !//////////////////////////////////////////////////////////////////////////// 
m = suavy(length(P), nprismsY);

fTV = 0.0;
 for i = 1:length(P) - 1
  fTV = fTV + sqrt(abs(P(i+1) - P(i))^2 + betaTV^2);
  fTV = fTV + sqrt(abs(P(m(i+1))-P(m(i)))^2 + betaTV^2);
 end

end