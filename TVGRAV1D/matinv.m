function [ ARRAY,SFLAG ] = matinv( ARRAY,NORDER)
% c	*************  Subrotina CHODEC *************************
% c
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
SFLAG=0;
DET=1; % 10
    for K=1:NORDER % 11
       AMAX=0;
       for I=K:NORDER
          for J=K:NORDER
            if (abs(AMAX)-abs(ARRAY(I,J))) <=0 % 24,24,30
              AMAX=ARRAY(I,J); % 24
              IK(K)=I;
              JK(K)=J;
             end
          end % 30
        end
       
      if (AMAX) == 0 % 41,32,41
       DET=0; % 32
       disp('AH,AH,AH, SE ESTREPOU: A MATRIZ E SINGULAR')
       SFLAG=1;
       return
      end

       I=IK(K); % 41
       if(I-K) >0 % 21,51,43
        for J=1:NORDER % 43
          SAVE1=ARRAY(K,J);
          ARRAY(K,J)=ARRAY(I,J);
          ARRAY(I,J)=-SAVE1; % 50
        end
       elseif (I-K) <0 %21
        continue;
       end
       
       J=JK(K); % 51
       
       if(J-K) > 0 % 21,61,53
        for I=1:NORDER % 53
         SAVE1=ARRAY(I,K);
         ARRAY(I,K)=ARRAY(I,J);
         ARRAY(I,J)=-SAVE1; % 60
        end
       elseif (J-K) < 0 % 21
        continue;
       end               
               
       for I=1:NORDER % 61
        if (I-K) ~= 0  % 63,70,63
         ARRAY(I,K)=-ARRAY(I,K)/AMAX; % 64
        end 
       end % 70

      for I=1:NORDER  % 71
          for J=1:NORDER
             if(I-K) ~=0 % 74,80,74
              if(J-K) ~=0 % 75,80,75  % 74
               ARRAY(I,J)=ARRAY(I,J)+ARRAY(I,K)*ARRAY(K,J); % 75 
              end
             end    
          end % 80
      end
      
       for J=1:NORDER % 81
        if(J-K) ~= 0 % 83,90,83
         ARRAY(K,J)=ARRAY(K,J)/AMAX; % 83
        end 
       end % 90

       ARRAY(K,K)=1./AMAX;
     end
      
DET=DET*AMAX; % 100
  for L=1:NORDER % 101
   K=NORDER-L+1;
   J=IK(K);
     if(J-K) >0 % 111,111,105
      for I=1:NORDER % 105
       SAVE1=ARRAY(I,K);
       ARRAY(I,K)=-ARRAY(I,J);
       ARRAY(I,J)=SAVE1; % 110
      end
    end
    I=JK(K); %111
   if(I-K) >0 % 130,130,113
     for J=1:NORDER % 113
      SAVE1=ARRAY(K,J);
      ARRAY(K,J)=-ARRAY(I,J);
      ARRAY(I,J)=SAVE1; % 120
     end
   end  
  end % 130
return;  % 140

end