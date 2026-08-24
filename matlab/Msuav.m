% Gerar matriz de vínculo por suavidade
function D = Msuav(Npv, Nph)
  D = [];
  for i=1:Nph
    for j=1:Npv-1
      M = zeros(Npv,Nph);
      M(j, i) = -1;
      M(j+1, i) = 1;
      D = [D; reshape(M, 1, Npv*Nph)];
    end
  end

  for i=1:Npv
    for j=1:Nph-1
      M = zeros(Npv,Nph);
      M(i, j) = -1;
      M(i, j+1) = 1;
      D = [D; reshape(M, 1, Npv*Nph)];
    end
  end
  D;
end