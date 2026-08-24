function f = funcGSS2D(m,x,grav_obs)
%FUNCGSS Summary of this function goes here
%   Detailed explanation goes here
% 
% global ab;
% global rhoa;
% global data;

     % forward TEST, output 
    dgrav1 =  Mod2DGravim(m);          % resistivitas hasil forward (calculated)
 lambda=0.05;

[J] = jacobian2DGrav(m,grav_obs,dgrav1);     
jac = inv(J'*J+lambda*eye(size(J'*J)));
dm = jac*J'*[(grav_obs)-(dgrav1)]';

mg = (m)+x*dm';
m=mg; 

% next model

   dgrav2 = Mod2DGravim(m);


f = ((norm(dgrav2-grav_obs))^2)/length(dgrav2);
end

