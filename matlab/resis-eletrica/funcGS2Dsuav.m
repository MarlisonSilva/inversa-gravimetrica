function f = funcGS2Dsuav(m,x,grav_cal,D)
%FUNCGSS Summary of this function goes here
%   Detailed explanation goes here
% 
% global ab;
% global rhoa;
% global data;

    

f =    (norm( log(Mod2DGrav(m))- (grav_cal))^2)/length(grav_cal) +x*norm(D*m')/length(m);




% next model

%   roac2 = Mod2DResist(m);


%f = ((norm(roac2-rhoa))^2)/length(roac2);
end
