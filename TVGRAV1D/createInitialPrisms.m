function [ P ] = createInitialPrisms( nobs,P,udist,anomObs,prismWidth,prismDensity )
%     !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% 	! createInitialPrisms
% 	!  
% 	! Brief:
% 	! 	For each prism, project its centre CP on the observations grid,
% 	! 	Define a circle C (observations grid) of radius r centered in CP
% 	! 	Get a weighted (1/DistanceFromCP**2) average of all the observations inside C
% 	! 	Use this value in the Bouguer approximation to find the prism thickness.
% 	!
% 	! Author: Williams Lima, williams_al@gmx.com
% 	!
% 	!//////////////////////////////////////////////////////////////////////////////////
   
ufact = 1e3; 
 if(udist == 0) % 0 - Distances in kilometers, 1 - distances in meters.
     ufact = 1.0;
 end				
dx = prismWidth;
rzero = dx/4.0;
r = rzero;
				
 for k = 1:size(P) 		 
  xcp = (k-1)*dx + dx/2;
  P(k) = 0.0;				 
   while 1
    m = 0;
    anom = 0.0;
    weight = 0.0;			
     for i = 1:nobs(1)						 
      d = abs( xcp - anomObs(i,1) );
      if((d < r) && (anomObs(i,2) ~= -9999.00))														 
	anom = anom + anomObs(i,2)/(d^2 + 0.0001);
	m = m + 1;
	weight = weight + 1.0/(d^2 + 0.0001);
      end
     end
				
  if(m ~= 0)
   anom = anom/weight;
   P(k) = ufact*anom/(41.9*prismDensity+0.0001);
   break;					 
  else
  % No points found. Increase the radius and try again
  r = r + rzero;
  end
  end
end

end