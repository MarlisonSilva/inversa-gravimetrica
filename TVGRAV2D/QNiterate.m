function [ chisqr, chisqA, YFIT, P, flambda ] = QNiterate( P, YFIT, anomObs, nobsX, nobsY, prismWidthX, prismWidthY, nprismsY, ...
                                                          prismDensity, dcType, dcBeta,  alpha, betaTV, flambda, iter,xmin, xmax, ymin, ymax, ...
                                                          nprisms_x, nprisms_y)
    
    nterms = length(P); 
    npts = nobsX*nobsY;   
    persistent  Hk;
    if isempty(Hk)
     Hk=eye(nterms,nterms); % Initialize Hessian matrix with identity.
    end
    
    persistent  bbeta;
    if isempty(bbeta)
     bbeta=zeros(nterms,1);
    end    

    idfp = 0;

%    nobs = length(anomObs(:,1));
%    lambda = flambda;
    %bbeta=0;
%    B=zeros(nterms,1);
    w=0;
        if(iter == 1)	
        %    ! Gradient of data at Pk.
            beta1 = zeros(nterms,1);		
            for i = 1:nobsX
                for j = 1:nobsY
                    w=j+(i-1)*nobsY; 
                    if(anomObs(w,3) ~=-9999.00)
                       DERIV = fderiv( i, j, P, anomObs, nobsX, nobsY, prismWidthX, prismWidthY, nprismsY, ...
                                     prismDensity);
                        beta1 = beta1 + (anomObs(w,3) - YFIT(i, j))*DERIV';
                    end 
                end
            end 
	     %   Gradient of total variation functional at Pk
            beta3 = gradFunctnTV(P, betaTV, nprismsY);
       
            %! Total gradient pk
            beta = -2*beta1 + alpha*beta3;
        else
            beta = bbeta;
        end    
        
    [chisq1, chisqA] = fchisq(anomObs, YFIT, P, alpha, betaTV, nprismsY, nobsX, nobsY);
    
    % Inversion of the modified curvature matrix to find new solution P.
    while 1
        
     for j = 1:nterms
      for k = 1:nterms
       ARRAY(j,k) = Hk(j,k) /sqrt(abs(Hk(j,j)*Hk(k,k)));
      end
       ARRAY(j,j) = ARRAY(j,j) + flambda;
     end
        beta=beta;
        %  matinv(ARRAY, nterms, nterms, DET, IK, JK, SFLAG)
%        ARRAY=chol(ARRAY);
        tic
        [ARRAY, SFlag]=CHODEC(ARRAY, nterms);
        toc
        for j=1:nterms
         for i=1:j
           ARRAY(i,j)=ARRAY(j,i);
         end
        end
            
        for j = 1:nterms
            B(j) = P(j);
            for k = 1:nterms
                B(j) = B(j)-beta(k)*ARRAY(j,k)/sqrt(abs(Hk(j,j)*Hk(k,k)));
                if(B(j) < 0)
                    B(j) = 0.01;
                end
            end
        end
  %      B=B';
    %    YFIT = functn(B,anomObs(:,1),xminPrisms,prismWidth,prismDensity,dcType,dcBeta,unit );
        YFIT = functn(B, anomObs, nobsX, nobsY, xmin, xmax, ymin, ymax, nprisms_x, nprisms_y, prismDensity);
        %chisqr = fchisq(anomObs(:,2),YFIT,B,alpha,betaTV,chisqA,tvValue);
      %  [chisqr,chisqA,tvValue] = fchisq( anomObs(:,2),YFIT,B,alpha,betaTV);
        [chisqr, chisqA] = fchisq(anomObs, YFIT, B, alpha, betaTV, nprismsY, nobsX, nobsY);
        
        if ((chisq1-chisqr) < 0)
            if(flambda > 1e8)
 %               flambda = -1.0;
                return;
            end
            flambda = 10.0*flambda;
            continue; %cycle
        else
            % Accepts solution, go ahead to the next iteration.
            % Gradient of data at Pk+1
            beta1 = zeros(nterms,1);
            for i = 1:nobsX
                for j = 1:nobsY
                    w=j+(i-1)*nobsY; 
                    if(anomObs(w,3) ~=-9999.00)
                       DERIV = fderiv( i, j, B', anomObs, nobsX, nobsY, prismWidthX, prismWidthY, nprismsY, ...
                                     prismDensity);
                        beta1 = beta1 + (anomObs(w,3) - YFIT(i, j))*DERIV';
                    end 
                end
            end 
            
            % Gradient of the total variation functional at Pk+1
            beta3 = gradFunctnTV(B', betaTV, nprismsY);
            
            % total gradient at pk+1
            bbeta = -2*beta1 + alpha*beta3;
            SK = B'-P;
            YK = bbeta-beta;
            
            if(idfp ~= 0)
                % Expression for DFP method start here.
                E1=0;
                for i = 1:nterms
                    E1 = E1+YK(I)*SK(I);
                end
                
                for i = 1:nterms
                    V1(I) = 0;
                    for j = 1:nterms
                        V1(i) = V1(i)+Hk(i,j)*SK(j);
                    end
                end
                
                for i = 1:nterms
                    V2(i) = YK(i)-V1(i);
                end
                
                E2=0;
                
                for i = 1:nterms
                    E2 = E2+V2(i)*SK(i);
                end
                
                for i = 1:nterms
                    for j = 1:nterms
                        XM1(i,j) = V2(i)*YK(j);
                    end
                end
                
                for i = 1:nterms
                    for j = 1:nterms
                        XM2(i,j) = YK(i)*YK(j);
                    end
                end
                
                for i = 1:nterms
                    for j = 1:nterms
                        XM3(i,j) =  (1./(E1+0.000001))*(XM1(i,j)+XM1(j,i))+ ...
                            (E2/(E1+0.000001)^2)*XM2(I,J);
                        Hk(i,j) = Hk(i,j)+XM3(i,j);
                        if(i.eq.j)
                            Hk(i,j) = hk(i,j);
                        end
                    end
                end
            else
                
                % Expression for BFGS starts here.
                
                for i = 1:nterms
                    for j = 1:nterms
                        XM2(i,j) = YK(i)*YK(j);
                    end
                end
                
                E1=0;
                
                for i = 1:nterms
                    E1 = E1+YK(i)*SK(i);
                end
                
                for i = 1:nterms
                    V1(i)=0;
                    for j = 1:nterms
                        V1(i) = V1(i)+Hk(i,j)*SK(j);
                    end
                end
                
                for i = 1:nterms
                    for j = 1:nterms
                        XM1(i,j) = V1(i)*V1(j);
                    end
                end
                
                E2=0;
                
                for i = 1:nterms
                    E2 = E2+SK(i)*V1(i);
                end
                
                for i = 1:nterms
                    for j = 1:nterms
                        XM3(i,j) = (1./(E1+0.000001))*XM2(i,j)-(1./(E2+0.000001))*XM1(i,j);
                        Hk(i,j) = Hk(i,j)+XM3(i,j);
                        if(i == j)
                            Hk(i,j) = Hk(i,j);
                        end
                    end
                end
            end
            
            % Update solution P
            P=B;
            flambda = flambda/10.0; 
            return
            
        end
    end
    
end