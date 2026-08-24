clear all;clc;
maxIter=100;
flambda=0.001;
usrTOL=0.0001;

noiseType=0;
noiseMean=0;
noiseSD=0;
noiseSeed=0;

alpha=1.25;
betaTV=0.001;

nprisms_x=42;
xmin=0;
xmax=52.5;
nprisms_y=26;
ymin=0;
ymax=32.5;
prismDensity=-0.30;
dcType=0;
dcBeta=0;
nobsX=21;
nobsY=13;

P=inv_2D( maxIter,flambda,usrTOL,noiseType,noiseMean,noiseSD,noiseSeed,alpha,betaTV,nprisms_x,nprisms_y,xmin, ...
                         xmax,ymin,ymax,prismDensity,dcType,dcBeta,nobsX,nobsY);
P=-1*P;       

clf('reset');
plotarmodelo2D(0,xmaxPrisms,P);