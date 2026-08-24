function [ pDeriv ] = gpDeriv(rho, xr, yr, hb, ht, xc, yc, xw2, yw2)
%     !C---------------------------------------------------------------
%     !C	DERIVADA ANALITICA DO PRISMA
%     !C---------------------------------------------------------------
%     !C	RO      = DENSIDADE
%     !C	XWIDTH2 = METADE DA LARGURA (EM X) DO PRISMA
%     !C	YWIDTH2 = METADE DA LARGURA (EM Y) DO PRISMA
%     !C	HTOP    = PROFUNDIDADE DO TOPO
%     !C	HBOTTON = PROFUNDIDADE DA BASE
%     !C	XCENTER = COORDENADA DO CENTRO DO PRISMA (EM X)
%     !C	YCENTER = COORDENADA DO CENTRO DO PRISMA (EM Y)
%     !C---------------------------------------------------------------
GAMMA = 6.67;
x1 = xc - xr - xw2;
x2 = xc - xr + xw2;
y1 = yc - yr - yw2;
y2 = yc - yr + yw2;
z1 = ht;
z2 = hb;

if(z2 == 0) 
 z2 = 1e3;
end
pDeriv = GAMMA*rho*( f_dkern(x1,y2,z2) + f_dkern(x2,y1,z2) -f_dkern(x2,y2,z2)  -f_dkern(x1,y1,z2) );

end