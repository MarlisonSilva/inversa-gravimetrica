function [ misfit ] = functn(P, XOBS, xminPrisms, prismWidth, prismDensity, dcType, dcBeta, unit)
%     !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%     ! functn
%     ! Compute the gravity field for a group of 2D justaposed prisms with fixed width.
%     ! density contrast constant or variable with depth.
%     !
%     ! Parameters:    
%     !   P               Input - Prisms thicknesses.
%     !   XOBS            Input - Observational points.
%     !   xminPrisms      Input - Start of the discretization.
%     !   prismWidth      Input - Prism width. The same for all prisms.
%     !   prismDensity    Input - Prisms contrast density. The same for all prisms.
%     !   dcType          Input -  Depth variation mode
%     !                           0 - Constant density contrast;
%     !                           1 - Hyperbolic variation with depth.    
%     !   dcBeta          Input - Density contrast decaying factor (used if dcType = 1).
%     !   unit            Input - 0 - Distances in kilometers. 1 - distances in meters.
%     ! Return:
%     !   Gravity field (mGal) at each point supplied in XOBS.
%     !
% 	!//////////////////////////////////////////////////////////////////////////////////////////////
  
% Parameters
 		
ztop = 0.0005;
misfit=zeros(length(XOBS),1);
for i = 1:length(P)
 xc = xminPrisms + (i-1)*prismWidth + prismWidth/2;			
 for j = 1:length(XOBS)
  switch dcType
    case 0
        % Constant density contrast
        misfit(j) = misfit(j) + grav2DPrism(prismDensity, xc, ztop, prismWidth, P(i), XOBS(j));                        
    case 1
        % Hyperbolic variation with depth
        misfit(j) = misfit(j) + hypGrav2DPrism(prismDensity, dcBeta, xc, ztop, prismWidth, P(i), XOBS(j));
  end 
 end 
end

        % Conversion from microGal to mGal if distances in meters.
if(unit == 1) 
  misfit = misfit/1000.0;
end  
end