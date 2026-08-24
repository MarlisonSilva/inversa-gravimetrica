function [ A, SFLAG ] = CHODEC(A,N)
% c	Esta subrotina calcula a inversa de ua matriz através da
% c	decomposição de Choleski. A matriz precisa ser simétrica 
% c	e definida-positiva.
% c
% c	ENTRADA:
% c
% c	A: matriz a ser invertida. Ela pode, opcionalmente ser 
% c	   definida como uma matriz triangular inferior. Se todos os 
% c	   elementos de A forem definidos, os elementos acima da 
% c	   diagonal principal serão preservados. Os elementos da 
% c	   diagonal e abaixo dela serão substituídos no output pelos
% c	   respectivos elementos da inversa de A.
% c	N: Ordem da matriz A.
% c	SFLAG: Indicador de erro:
% c		   SFLAG=0: A inversa foi calculada corretamente
% c		   SFLAG=1: O cálculo foi interrrompido e o controle
% c		      transferido para o programa principal.
% c
% c	SAÍDA:
% c
% c	A: Matriz N x N contendo a inversa de A nos elementos diagonais
% c	   e abaixo dela. Os elementos acimas da diagonal não são 
% c	   modificados por esta subrotina.
% c
% c
%A=zeros(N,N);
IFLAG=0;
SFLAG=0;

    if (A(1,1) < 0) %GO TO 20
     IFLAG=1;
     for J=1:N
      for I=J:N
       A(I,J)=-A(I,J); % 10
      end
     end
    end
   A(1,1)=sqrt(A(1,1)); % 20
   
 for IR=2:N
  for II=1:IR
   SUM=0;
  if(II ~= IR)
   if(II ~= 1)
     for J=1:II-1
      SUM=SUM+A(IR,J)*A(II,J);
     end 
   end  
   A(IR,II)=(A(IR,II)-SUM)/A(II,II);
   continue;% break; % GO TO 90
  end 

  for J=1:IR-1
   SUM=SUM+A(IR,J)*A(IR,J);
  end
  A(IR,II)=A(IR,II)-SUM;
 % IF(A(IR,II).LT.0.) GO TO 70
  A(IR,II)=sqrt(A(IR,II));
%	break; % GO TO 90
%70	print 80
%80	FORMAT(10X,'SUBROTINA INTERROMPIDA ',/,10X,'A DECOMPOSICAO DE
%	1CHOLESKI',
%    *  ' SO PODE SER APLICADA A MATRIZES POSITIVO-DEFINIDAS.',
%    * /,10X,'TENTE A REDUCAO DE GAUSS-JORDAN.')
%	SFLAG=1;
%	return;
  end % 90
 end

 for J=1:N
  for I=J:N
  if(I == J)
   A(I,I)=1./A(I,I); 
 else
   SUM=0;
    for K=J:I-1
     SUM=SUM+A(I,K)*A(K,J);
    end
    A(I,J)=-SUM/A(I,I);
   end
  end
 end
 %	GO TO 120
 
%110	IF(A(I,I).EQ.0.) GO TO 130
%	A(I,I)=1./A(I,I);
%120	CONTINUE
%	GO TO 145
%130	print 140
%140	FORMAT (10X,'CHODEC: A MATRIZ E SINGULAR.')
%	SFLAG=1.
%return
 
 for J=1:N 
  for I=J:N
   SUM=0;
    for K=I:N
     SUM=SUM+A(K,I)*A(K,J);
    end
   A(I,J)=SUM;
  end 
 end

 if (IFLAG == 0) 
  return;
 end 
 
 for J=1:N
  for I=J:N
	A(I,J)=-A(I,J);
  end
 end

 return;
end