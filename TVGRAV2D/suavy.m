function [ m ] = suavy(n,nx)
k=0;
 for j=1:nx
  for i=j:nx:n
   k=k+1;
   m(k)=i;
  end
 end
return;
end