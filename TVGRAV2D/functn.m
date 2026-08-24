function [ Grav ] = functn(P, anomObs, nobsX, nobsY, xmin, xmax, ymin, ymax, nprismsX, nprismsY, prismDensity)
%     ! Compute the gravity field for a group of juxtaposed 3D prisms    
%     !
%     ! Parameters:    
%     !   P               Input - Prisms thicknesses.
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
%     ! Return:
%     !   Gravity field (mGal) at each point on the observations grid
%     !
%     !
% 	!//////////////////////////////////////////////////////////////////////////////////////////////
htop = 0.0001;

% nobsX = size(anomObs(:,1,1))
% nobsY = size(anomObs(1,:,1))

dx = (xmax-xmin)/nprismsX;
dy = (ymax-ymin)/nprismsY;

%dx = prismWidthX;
%dy = prismWidthY;

Grav=zeros(nobsX,nobsY);		
w=0;        
for i = 1:nobsX
    for j = 1:nobsY
        w=j+(i-1)*nobsY; 
        xob = anomObs(w,1);
        yob = anomObs(w,2);
        valor=anomObs(w,3);
        if(anomObs(w,3) ~= -9999.00)
            for k = 1: length(P)
                row = fix((k-1)/nprismsY + 1);
                col = fix(k - (row - 1)*nprismsY);
                hbase = P(k) + htop;
                xcr = ((row-1)*dx + row*dx)/2;								
                ycr = ((col-1)*dy + col*dy)/2;
                Grav(i,j) = Grav(i,j) + GPRISM(xob, yob, dx/2, dy/2, htop, ...
                                        hbase, xcr, ycr, prismDensity);
            end
        else
            Grav(i,j) = -9999.00;
        end    
    end
end
end