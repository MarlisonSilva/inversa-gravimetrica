clear all
clc
close all

load('gobsGI.mat');
load('modeloGI.mat');
%load Dgrav  % Dado observado
%load Xini % modelo inicial
%


%modGrav1d2( P, XOBS, iobs, xminPrisms, prismWidth, prismDensity )            
Dgrav = gobsGI ;

modeloinicial=Dgrav*100;

m=modeloinicial;
lr = length(m); 
kr = 1e-4;                  % tolerancia
j = 1;                           % initial iteration
iteration(j) = 1;
itermax = 5; % 100

grav_cal= Mod2DGravim(m);


rms_err = norm(grav_cal-Dgrav)/sqrt(length(Dgrav));
rms_error(j) = rms_err;



while(rms_err > kr)
    
    % find lamda value using golden section search
    % m,a,b,grav_cal
    [lamda]=0.1;
    % [lamda] = gss2D_lm(m,0.001,20,grav_cal);  
   % [delta] = gss2D_lm(m,0.001,1,grav_cal);
   delta=1;
   %  disp('lamda'); disp(lamda);
   %  lamda_plot(j) = lamda;

    % === Find optimum lamda === %
    % Jacobian matrix 
    [J] = jacobian2DGrav(m,Dgrav,grav_cal);    
    
    % Levenberg-marquardt algorithm
    jac = inv(J'*J+lamda*eye(size(J'*J)));
    dm = jac*J'*[Dgrav-(grav_cal)]';
    mg = (m)+delta*dm';
    m=mg;
     
     
    % find next model
   
       grav_cal =Mod2DGravim(m);
   
     
    rms_err = norm(grav_cal-Dgrav)/sqrt(length(Dgrav));
   Mhist(:,:,j)=-m;
    j = j + 1;
    rms_error(j) = rms_err; 
    iteration(j) = iteration(j-1)+1; 
     if (iteration(j) > itermax)
         break
     end
    % plot_iterasi(m,rhoa_cal,j);
 
end


  
  %%
  xobs=1:60;
  plot(xobs, -modeloGI, 'LineWidth',3);
  hold on
   plot(xobs, m, 'LineWidth',3);
  




