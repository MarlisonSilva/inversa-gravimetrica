% inversao 2d resistividade % 
% Programador Francisco Marcio Barboza %
clear all
clc
close all; 

loadlibrary('rho2D.dll','rho2D.h');

%%


clc; 

load X % modelo verdadeiro

 roav=Mod2DResist(X); %dado calculado sem ruido
 
 roa=roav+0.025*randn(size(roav)).*roav; % dado observado com ruido
 
 Xinic=mean(X)*ones(size(X));% modelo inicial tipo homogeneo
 %%
m = Xinic;
 
lr = length(m); 
kr = 10e-10; 
iteration = 1; 
maxiteration = 6; 
dfit = 1;


while iteration<maxiteration

roa1=Mod2DResist(m);

e1 = [log(roa)-log(roa1)];
dd = e1; 
misfit1 = e1' *e1;
if misfit1<kr
break
end
[A] = jacobian2D(m,roa,roa1);

 say = 1;
 k = 0;

dmg=inv(A'*A+0.01*eye(size(A'*A)))*A'*dd';% depois implemento suavidade
% dmg=pinv(A'*A)*A'*dd';
mg = exp(log(m)+dmg);
roa4=Mod2DResist(mg);

e2 = [log(roa)-log(roa4)];
misfit2 = e2'*e2;
if misfit2>misfit1
('Beta control');
say = say+1; 
k = k+1;

else
m = mg; 
dfit = (misfit1-misfit2)/misfit1;
iteration = iteration+1; 
a = iteration;
if dfit<kr
iteration = maxiteration; 
say = say+1;
end
end
end

observed = roa; 
calculated = roa4;
m;

%%
Mest=reshape(m,18,185);

imagesc(Mest) % imagem do modelo so para ter noção


