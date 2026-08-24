function [ deriv ] = fderiv(i, j, P, anomObs, nobsX, nobsY, prismWidthX, prismWidthY, nprismsY, ...
                        prismDensity)
%     !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%     !   Computes the partial derivatives of the fitting functional with respect to Pj at Xij.
%     !   Where Pj is th jth element of solution vector P and Xij is the position of the
%     !   (i,j)th observation.
%     !
%     ! Parameters:  
%     !   i,j	            Input - Observation index.
%     !   P               Input - Solution vector (prisms thicknesses)
%     !   anomObs         Input - Observations
%     !   prismWidthX     Input - Horizontal prisms width along X (N-S)
%     !   prismWidthY     Input - Horizontal prisms width along Y (E-W)
%     !   nprismsY        Input - Number of prisms along Y (E-W)
%     !   prismDensity    Input - Prisms density contrast. The same for all prisms.
%     !   dcType          Input - 
%     !                           0 - Constant density contrast
%     !                           1 - Hyperbolic variation with depth
%     !                           2 - Parabolic variation with depth
%     !   dcBeta          Input - Density contrast decaying factor    
%     !   DERIV           Output - Derivatives for each Pj
%     !
%     !////////////////////////////////////////////////////////////////////////////////////////

deriv = 0.0; 
w=0;
 for k = 1:size(P)            
 row = fix((k-1)/nprismsY + 1);
 col = fix(k - (row - 1)*nprismsY);
        
 xcr = ((row-1)*prismWidthX + row*prismWidthX)/2;								
 ycr = ((col-1)*prismWidthY + col*prismWidthY)/2;
 w=j+(i-1)*nobsY;
 deriv(k) = gpDeriv( prismDensity, anomObs(w,1),anomObs(w,2), 0.0001+P(k),0.0001, ...
            xcr,ycr,prismWidthX/2,prismWidthY/2 );

 end

end