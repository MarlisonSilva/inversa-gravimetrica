function [ P ] = createInitialPrisms( P,prismDensity,anomObs,nobsX,nobsY,nprisms_x, nprisms_y, xmax, xmin, ymax, ymin )
%   !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%   ! createInitialPrisms
%   !  
%   ! Brief:
%   !   For each prism, project its centre CP on the observations grid,
%   !   Define a circle C (observations grid) of radius r centered in CP
%   !   Get a weighted (1/DistanceFromCP**2) average of all the observations inside C
%   !   Use this value in the Bougher approximation to find the prism thickness.
%   !
%   !//////////////////////////////////////////////////////////////////////////////////

EPS = 1E-4;
    
dx = (xmax-xmin)/nprisms_x;
dy = (ymax-ymin)/nprisms_y;
rzero = dx/4;
r = rzero;

for k = 1:length(P)
       
 row = fix((k-1)/nprisms_y + 1);
 col = fix(k - (row - 1)*nprisms_y);
      
 xcp = ((row-1)*dx + row*dx)/2;
 ycp = ((col-1)*dy + col*dy)/2;       
       
 P(k) = 0.0;         
 while 1
  m = 0;
  campo = 0;
  peso = 0;
  w=0;
   for i = 1:nobsX
     for j = 1:nobsY
        w=j+(i-1)*nobsY; % w=i+(j-1)*nobsX; 
        d = sqrt((xcp - anomObs(w,1))^2 + (ycp - anomObs(w,2))^2);
        if((d < r) && (anomObs(w,3) ~= -9999.00))				                       
          campo = campo + anomObs(w,3)/(d^2+EPS);
          m = m + 1;
          peso = peso + 1./(d^2+EPS);
        end
     end
   end

          if(m == 0) 
           % ! No points found. Increase the radius and try again'
            r = r + rzero;
          else
	    campo = campo/peso;
	    P(k) = campo/(41.9*prismDensity+EPS);
           %  if k==13
          %    r=r;   
         %    end  
            break; %exit;
          end
  
 end
end
end