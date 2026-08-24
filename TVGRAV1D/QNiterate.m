function [ chisqr, chisqA, tvValue,YFIT, P ] = QNiterate( nprisms,P,anomObs,YFIT,xminPrisms,prismWidth,prismDensity,dcType,dcBeta, ...
        alpha,betaTV,flambda,iter,unit)
    
    nterms = length(P);    
    persistent  Hk;
    if isempty(Hk)
     Hk=eye(nprisms,nprisms); % Initialize Hessian matrix with identity.
    end
    
    persistent  bbeta;
    if isempty(bbeta)
     bbeta=zeros(nterms,1);
    end    

    idfp = 0;

    nobs = length(anomObs(:,1));
    lambda = flambda;
    %bbeta=0;
    B=zeros(nterms,1);
    
    if(iter == 1)
        % Gradient of data at Pk.
        beta1 = 0;
        for i = 1: nobs(1)
            DERIV=fderiv(i,P,anomObs(:,1),xminPrisms,prismWidth,prismDensity,dcType,dcBeta,unit);
            beta1 = beta1 + ( anomObs(i,2) - YFIT(i) )*DERIV;
        end
        % Gradient of total variation functional at Pk
        beta3 = gradFunctnTV(P, betaTV);
        % Total gradient pk
        beta = -2*beta1 + alpha*beta3';
     else
        beta = bbeta;
    end
    
    [chisq1,chisqA,tvValue] = fchisq( anomObs(:,2),YFIT,P,alpha,betaTV);
    
    % Inversion of the modified curvature matrix to find new solution P.
    while 1
        
        for j = 1:nterms
            for k = 1:nterms
                ARRAY(j,k) = Hk(j,k) /sqrt(abs(Hk(j,j)*Hk(k,k)));
            end
            ARRAY(j,j) = ARRAY(j,j) + lambda;
        end
        
        %  matinv(ARRAY, nterms, nterms, DET, IK, JK, SFLAG)
%        ARRAY=chol(ARRAY);
        [ARRAY, SFlag]=matinv(ARRAY, nterms);
        for j = 1:nterms
            B(j) = P(j);
            for k = 1:nterms
                B(j) = B(j)-beta(k)*ARRAY(j,k)/sqrt(abs(Hk(j,j)*Hk(k,k)));
                if(B(j) < 0)
                    B(j) = 0.01;
                end
            end
        end
        
        YFIT = functn(B,anomObs(:,1),xminPrisms,prismWidth,prismDensity,dcType,dcBeta,unit );
        %chisqr = fchisq(anomObs(:,2),YFIT,B,alpha,betaTV,chisqA,tvValue);
        [chisqr,chisqA,tvValue] = fchisq( anomObs(:,2),YFIT,B,alpha,betaTV);
        
        if ((chisq1-chisqr) < 0)
            if(lambda > 1e8)
                flambda = -1.0;
                return
            end
            lambda = 10.0*lambda;
            continue; %cycle
        else
            % Accepts solution, go ahead to the next iteration.
            % Gradient of data at Pk+1
            beta1 = 0;
            for i = 1:nobs
                DERIV=fderiv(i,B,anomObs(:,1),xminPrisms,prismWidth,prismDensity,dcType,dcBeta,unit);
                beta1 = beta1 + (anomObs(i,2) - YFIT(i))*DERIV;
            end
            % Gradient of the total variation functional at Pk+1
            bbeta3 = gradFunctnTV(B, betaTV);
            % total gradient at pk+1
            bbeta = -2*beta1' + alpha*bbeta3;
            SK = B-P;
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
            return
            
        end
    end
    
end