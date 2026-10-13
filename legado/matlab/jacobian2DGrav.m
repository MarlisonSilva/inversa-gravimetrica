function[A] = jacobian2DGrav(m,Dgrav,grav_cal)

lr=length(m);


par = 0.1;
r2 = m;
for i2 = 1:lr
r2(i2) = (m(i2)*par)+m(i2);


Dgrav2= Mod2DGravim( r2);



A1(:,i2) = ((Dgrav2-grav_cal)/(m(i2)*par))*m(i2)./Dgrav;
r2 = m;
end
A = A1;
return
end