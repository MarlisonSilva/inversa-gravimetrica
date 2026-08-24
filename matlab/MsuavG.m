
function D = MsuavG(Np)
% Matriz n-1 por n de suavidade 
%para o problema relevo do embasamento 2D
% np= numero de paramentros
% Matriz 
D = [];
    for i=1:Np-1
      M = zeros(1,Np);
      M(i) = -1;
      M(i+1) = 1;
      D(i,:)=M;
    end


  D;
end