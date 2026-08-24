clear all
clc
close all; 

loadlibrary('rho2D.dll','rho2D.h');

%unloadlibrary rho2D;

%%
%global rhoa;

%%
clear
clc
load Xini
load roa

rhoa=roa;

%%
Nph=85;
Npv=18;

D = Msuav(Npv, Nph);

%%

m=Xinic/5;
lr = length(m); 
kr = 10e-11;                  % tolerancia
j = 1;                           % initial iteration
iteration(j) = 1;
itermax = 7; % 100

rhoa_cal= Mod2DResist(m);

rms_err = norm(rhoa_cal-rhoa)/sqrt(length(rhoa));
rms_error(j) = rms_err;



while(rms_err > kr)
    
    % find lamda value using golden section search
     [lamda] = gss2D_lmsuav(m,0.001,20,rhoa,D);  
     disp('lamda'); disp(lamda);
     lamda_plot(j) = lamda;

    % === Find optimum lamda === %
    % Jacobian matrix (using input: X, rhoa_cal)
    [J] = jacobian2D(m,rhoa,rhoa_cal);    
    
    % Levenberg-marquardt algorithm
    jac = inv(J'*J+lamda*(D'*D));
    dm = jac*J'*[log(rhoa)-log(rhoa_cal)]';
    mg = exp(log(m)+dm');
    m=mg;
     
     
    % find next model
   
       rhoa_cal =Mod2DResist(m);
   
     
    rms_err = norm(rhoa_cal-rhoa)/sqrt(length(rhoa));
   Mhist(:,:,j)=m;
    j = j + 1;
    rms_error(j) = rms_err; 
    iteration(j) = iteration(j-1)+1; 
     if (iteration(j) > itermax)
         break
     end
    % plot_iterasi(m,rhoa_cal,j);
    
end





%%

load('X.mat')
mverd= reshape(X,18,185);

mverd=mverd(:,51:135);
Mest=reshape(m,18,85);
figure(1)
imagesc(Mest) 
figure(2)
imagesc(mverd) 

