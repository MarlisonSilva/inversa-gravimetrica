function [ P ] = inv_2D( maxIter,flambda,usrTOL,noiseType,noiseMean,noiseSD,noiseSeed,alpha,betaTV,nprisms_x,nprisms_y,xmin, ...
                         xmax,ymin,ymax,prismDensity,dcType,dcBeta,nobsX,nobsY)

load anomObs.mat;                     
%prismWidth = (xmaxPrisms - xminPrisms)/nprisms;

P=zeros(nprisms_x*nprisms_y,1);
nobs = size(anomObs);
YFIT=zeros(nobs(1),1);
NOISEV=zeros(nobs(1),1);


%NOISEV=noise(noiseSeed, nobs(1), noiseMean, noiseSD, noiseType);

mu=0;
sigma=1/80; % Adiciona ruído gaussiano em torno de 5%
%NOISEV=sigma*randn(nobs(1),1)+mu;
%anomObs(:,3,:) = anomObs(:,3) + anomObs(:,3).*NOISEV;
P=createInitialPrisms(P,prismDensity,anomObs,nobsX,nobsY,nprisms_x, nprisms_y, xmax, xmin, ymax, ymin);

YFIT = functn(P, anomObs, nobsX, nobsY, xmin, xmax, ymin, ymax, nprisms_x, nprisms_y, prismDensity);

iter=1;
convCont =0;
dx = (xmax-xmin)/nprisms_x;
dy = (ymax-ymin)/nprisms_y;

while 1
 [chisqr, chisqA, YFIT, P, flambda] = QNiterate( P, YFIT, anomObs, nobsX, nobsY, dx, dy, nprisms_y, ...
                                                          prismDensity, dcType, dcBeta,  alpha, betaTV, flambda, iter, xmin, xmax, ymin, ymax, ...
                                                          nprisms_x, nprisms_y);

 % Termination criteria

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

disp([iter, chisqr, chisqA]);
chisq1 = chisqr; 

 if(iter == maxIter)
  disp('Maximum number of iterations reached!');
  break;
 end
 
iter=iter+1;
rnit = iter;

end
 
 
end