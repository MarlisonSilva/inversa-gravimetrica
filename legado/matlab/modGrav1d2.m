function functn= modGrav1d2( P, XOBS, iobs, xminPrisms, prismWidth, prismDensity )            
%        implicit none
    
%         ! Parameters
%         real, dimension(:), intent(in) :: P, XOBS
% 		real, intent(in)	:: xminPrisms, prismWidth, prismDensity, dcBeta
% 		integer, intent(in)	:: dcType, unit
% 		
%         ! Local variables
%         integer :: i,j             	    
% 		real	:: xc, ztop = 0.0005
%         real, dimension(5)  :: X, Z
%         real    :: tmp
% 
% 		! Return
% 		real, dimension(size(XOBS)) :: functn
% 
% 		!\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%         ! Start of execution block
%         !////////////////////////////
      
		functn=0;
        
         ztop = 0.0005;
		for i = 1: length(P)
		    xc = xminPrisms + (i-1)*prismWidth + prismWidth/2	;		
			%for j = 1: length(XOBS)
			    %select case (dcType)
%                 switch dcType
%                 case {'0'}
                       % ! Constant density contrast
                         functn = functn + grav2DPrism2(prismDensity, xc, ztop, prismWidth, P(i), XOBS(iobs))  ;     % grav2DPrism(den, x0, z0,w, h, xobs)                  
%                     case {'1'}
%                        % ! Hyperbolic variation with depth
%                         functn(:,j+1) = functn(:,j) + hypGrav2DPrism(prismDensity, dcBeta, xc, ztop, prismWidth, P(i), XOBS(j));
%                end 
            %end 
       %  FF(i)=functn;   
		end 
%F= functn;
       % ! Conversion from microGal to mGal if distances in meters.
	%	if(unit==1) 
            functn = functn/1000.0;
    %    end

    end