function [ DERIV ] = fderiv( i,P,XOBS,xminPrisms,prismWidth,prismDensity,dcType,dcBeta,unit)
%     !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%     ! fderiv
%     !   Computes the partial derivatives of the fitting functional with respect to Pj at Xi.
%     !   Where Pj is the jth element of solution vector P and Xi is the position of the
%     !   ith observation.
%     !
%     ! Parameters:  
%     !   i	            Input - Observation index.
%     !   P               Input - Solution vector (prisms thicknesses)
%     !   XOBS            Input - Observational points
%     !   xminPrisms      Input - Start of the discretization 
%     !   prismWidth      Input - Prism width. The same for all prisms.
%     !   prismDensity    Input - Prisms contrast density. The same for all prisms.
%     !   unit            Input - 0 - Distances in kilometers. 1 - distances in meters.
%     !   DERIV           Output - Derivatives for each Pj
%     !
%     !////////////////////////////////////////////////////////////////////////////////////////
    
switch dcType
case 0
    % Constant density contrast
    for k = 1:size(P)
        xc = xminPrisms + (k-1)*prismWidth + prismWidth/2;
        DERIV(k) = dgrav2DPrism(prismDensity, xc, 0.0005, prismWidth, P(k), XOBS(i), unit);
    end
case 1
    % Hyperbolic density contrast
    for k = 1:size(P)
        xc = xminPrisms + (k-1)*prismWidth + prismWidth/2;
        DERIV(k) = dhypGrav2DPrism(prismDensity, dcBeta, xc, 0.0005, prismWidth, P(k), XOBS(i));
    end
end  


end