function f = funcGS2D(m,x,Dgrav)
%FUNCGSS Summary of this function goes here
%   Detailed explanation goes here
% 
% global ab;
% global rhoa;
% global data;

     % forward TEST, output 
    %roac1 =  Mod2DResist(m);          % resistivitas hasil forward (calculated)


%[J] = jacobian2D(m,rhoa,roac1);     
%f =    inv(J'*J+x*eye(size(J'*J)));

f =    (norm((Mod2DGravim(m))- (Dgrav))^2)/length(Dgrav) +x*norm(m)/length(m);





end

