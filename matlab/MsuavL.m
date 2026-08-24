function D = MsuavL(Npv, Nph)
  D = [];
  for i=1:Nph
    for j=1:Npv-1
      M = zeros(Npv,Nph);
      M(j, i) = -1;
      M(j+1, i) = 1;
      D = [D; reshape(M, 1, Npv*Nph)];
    end
  end
  D;
end