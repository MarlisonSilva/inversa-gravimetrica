function [ A ] = plotarmodelo2D( xminPrisms,xmaxPrisms,P )

nprisms=length(P);
prismWidth = (xmaxPrisms - xminPrisms)/nprisms;

x1=xminPrisms;
x2=prismWidth;
y1=0;

for i=1:nprisms
 y2=P(i);

 x = [x1, x2, x2, x1, x1];
 y = [y1, y1, y2, y2, y1];
 plot(x, y, 'b-');
 
 x1=x2;
 x2=x1+prismWidth;
 
 hold on;
end
 xlim([ xminPrisms, xmaxPrisms]);
 ylim([min(P)+0.5*min(P), 0]);

% x1=0.5;
% x2=1;
% 
% y1=0;
% y2=-1;
% 
% x = [x1, x2, x2, x1, x1];
% y = [y1, y1, y2, y2, y1];
% plot(x, y, 'b-');

end