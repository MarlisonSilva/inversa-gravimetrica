function [ P ] = inv_1D( maxIter,flambda,usrTOL,noiseType,noiseMean,noiseSD,noiseSeed,alpha,betaTV,nprisms,xminPrisms, ...
                         xmaxPrisms,prismDensity,dcType,dcBeta,udist)

load anomObs.mat;                     
prismWidth = (xmaxPrisms - xminPrisms)/nprisms;

P=zeros(nprisms,1);
nobs = size(anomObs);
YFIT=zeros(nobs(1),1);
NOISEV=zeros(nobs(1),1);


%NOISEV=noise(noiseSeed, nobs(1), noiseMean, noiseSD, noiseType);

mu=0;
sigma=1/80; % Adiciona ruído gaussiano em torno de 5%
NOISEV=sigma*randn(nobs(1),1)+mu;
anomObs(:,2) = anomObs(:,2) + anomObs(:,2).*NOISEV;
P=createInitialPrisms(nobs,P,udist,anomObs,prismWidth,prismDensity);

YFIT = functn(P,anomObs(:,1),xminPrisms,prismWidth,prismDensity,dcType,dcBeta,udist);

iter=0;
convCont =0;

while 1
 [ chisqr, chisqA, tvValue,YFIT, P  ] = QNiterate( nprisms,P,anomObs,YFIT,xminPrisms,prismWidth,prismDensity,dcType,dcBeta, ...
                            alpha,betaTV,flambda,iter+1,udist);
 iter=iter+1;
 rnit = iter;

 disp([iter, chisqr, chisqA]);

 % Termination criteria
 if(flambda < 0) 
   return
 end

 if(iter > 1) 		
  conv = (abs(chisq1-chisqr))/chisq1;			            
  if(conv <= usrTOL) 					
   convCont = convCont + 1;
   if(convCont >= 5)					
    break %return
   end
  else
   convCont = 0; 
  end					
 end


 if(iter == maxIter)
  disp('Maximum number of iterations reached!');
  break;
 end

 
chisq1 = chisqr;
end
 
 
end