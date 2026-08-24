function [ prism ] = GPRISM(XR,YR,XW2,YW2,HT,HB,XCR,YCR,RHO)
GAMA = 6.67;
XC   = XCR - XR;
YC   = YCR - YR;
X1   = XC - XW2;
X2   = XC + XW2;
Y1   = YC - YW2;
Y2   = YC + YW2;
Z1   = HT;
Z2   = HB;
prism=GAMA*RHO*(FKERN(X1,Y1,Z1)+FKERN(X1,Y2,Z2) + ...
         FKERN(X2,Y1,Z2)+FKERN(X2,Y2,Z1)-FKERN(X2,Y2,Z2)- ...
         FKERN(X2,Y1,Z1)-FKERN(X1,Y2,Z1)-FKERN(X1,Y1,Z2));
    
return;
end