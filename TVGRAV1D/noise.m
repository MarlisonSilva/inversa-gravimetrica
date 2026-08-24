function [ r ] = noise( iy,n,xm,sd,IDIST)
% 	!	Esta subrotina gera uma sequência de N números pseudoaleatórios
% 	!	com distribuições uniforme, Gaussiana, Laplace ou Cauchy.
% 	!
% 	!		INPUT:
% 	!
% 	!		IY:    semente para inicialização da sequência aleatória
% 	!		N:     número de pontos da sequência;
% 	!		XM:    média da distribuição;		
% 	!		SD:    desvio padrão no caso de  distribuição Gaussiana e 
% 	!			   de Laplace;
% 	!			   semi-largura no caso da distribuição de Cauchy;
% 	!			   metade do intervalo em que a fdp é diferente de zero
% 	!			   no caso da distribuição uniforme;
% 	!		IDIST: parâmetro de escolha da distribuição:
% 	!			    =1: distribuição uniforme
% 	!			    =2: distribuição Gaussiana
% 	!			    =3: distribuição de Laplace
% 	!			    =4: distribuição de Cauchy
% 	!
% 	!		OUTPUT:
% 	!
% 	!		R:     vetor de dimensão N contendo a sequência aleatória
% 	!
% 	! Author: Prof. Joao Batista
% 	!
r=zeros(n,1);

for i=1:n
 switch IDIST	
  case 1
   r(i)=urandu(iy,xm,sd);
  case 2
   r(i)=urandn(iy,xm,sd);
  case 3
   r(i)=urandl(iy,xm,sd);
  case 4
   r(i)=urandc(iy,xm,sd);
 end
end

end

function [X] = urandu(IY,XM,SD)
% 	!	Esta subrotina gera uma realização X de uma variável 
% 	!	pseudo-aleatória Gaussiana com média XM e desvio padrão SD. 
% 	!	IY é a semente (inteiro) que inicializa uma sequência particular, 
% 	!	reproduzível, de realizações. A subrotina se encarrega de gerar 
% 	!	automaticamente outros valores para a semente, cada vez que 
% 	!	ela é chamada.
% 	!
% 	!	Exemplo de uso:
% 	!
% 	!	n=100
% 	!	IY=9
% 	!	do  i=1,n
% 	!	call urandn(IY,0.,1.,r)
% 	!	y(i)=y(i)+r
% 	!	end do
% 	!
% 	!	Função necessária: URAND(IY)
% 	!
% 	!	Author: Prof. Joao Batista
% 	! 
S=10;
 while S>1
  U1=urand(IY);
  U2=urand(IY);
  V1=2*U1-1;
  V2=2*U2-1;
  S=V1*V1+V2*V2;
%	IF(S.GT.1.) GO TO 1
 end
 P=SQRT(-2.*ALOG(S)/S)*SD;
 X=V1*P+XM;
end  


function [X] = urand(IY)
% 	!	Esta função gera uma realização de uma variável pseudo-aleatória
% 	!	com distribuição uniforme no intervalo [0,1].
% 	!	IY é a semente (inteiro) que inicializa uma sequência particular, 
% 	!	reproduzível, de realizações. A função se encarrega de gerar 
% 	!	automaticamente outros valores para a semente, cada vez que 
% 	!	ela é chamada.
% 	!
% 	!	Exemplo de uso:
% 	!
% 	!	n=100
% 	!	iy=3
% 	!	do  i=1,n
% 	!	y(i)=y(i)+urand(iy)
% 	!	end do
% 	!
% 	!	Author: Prof. Joao Batista
% 	!   

%M2=0;
ITWO=2;
%	IF(M2.NE.0) GO TO 20
M=1;
M2=M;
M=ITWO*M2;
%IF(M > M2) GO TO 10
HALFM=M2;
IA=8*int64(HALFM*atan(1)/8)+5;
IC=2*int64(HALFM*(0.5-sqrt(3)/6))+1;
MIC=(M2-IC)+M2;
S=0.5/HALFM;
IY=IY*IA;

 if(IY > MIC) 
  IY=(IY-M2)-M2;
 end
 IY=IY+IC;
 if(IY/2 > M2) 
  IY=(IY-M2)-M2;
 end
 if(IY < 0) 
  IY=(IY+M2)+M2;
 end 
 X=float(IY)*S;
end   

function [X] = urandn(IY,XM,SD)
% 	!	Esta subrotina gera uma realização X de uma variável 
% 	!	pseudo-aleatória Gaussiana com média XM e desvio padrão SD. 
% 	!	IY é a semente (inteiro) que inicializa uma sequência particular, 
% 	!	reproduzível, de realizações. A subrotina se encarrega de gerar 
% 	!	automaticamente outros valores para a semente, cada vez que 
% 	!	ela é chamada.
% 	!
% 	!	Exemplo de uso:
% 	!
% 	!	n=100
% 	!	IY=9
% 	!	do  i=1,n
% 	!	call urandn(IY,0.,1.,r)
% 	!	y(i)=y(i)+r
% 	!	end do
% 	!
% 	!	Função necessária: URAND(IY)
% 	!
% 	!	Author: Prof. Joao Batista
% 	! 
S=1;
 while S>1
 U1=urand(IY);
 U2=urand(IY);
 V1=2.*U1-1;
 V2=2.*U2-1;
 S=V1*V1+V2*V2;
%	IF(S.GT.1.) GO TO 1
 end 
P=sqrt(-2.*log(S)/S)*SD;
X=V1*P+XM;
end  

function [X] = urandl(IY,VM,DL)
% 	!
% 	!
% 	!	SUBROTINA URANDL GERA NUMEROS ALEATORIOS COM DISTRIBUICAO
% 	!   DE LAPLACE COM MEDIA=VM E PARAMETRO DE DISPERSAO = DL
% 	!
% 	!	Author: Prof. Joao Batista
% 	!
W=urand(IY);
SL=DL;
 if(W <= 0.5) 
  X=SL*ALOG(2*W)+VM;
 else
  SL=-DL;
  W=1-W;
  X=SL*ALOG(2*W)+VM;  
 end

end    

function [X] = urandc(IY,VM,SL)
% 	!
% 	!	SUBROTINA URANDC GERA NUMEROS ALEATORIOS COM
% 	!  DISTRIBUICAO DE CAUCHY COM VALOR MAIS FREQUENTE =
% 	!  VM E SEMI-LARGURA=SL
% 	!
% 	! Author: Prof. Joao Batista
% 	! 
U=urand(IY);
PI=4*atan(1);
ALPHA=PI*(2*U-1)/2;
SALF=sin(ALPHA);
CALF=cos(ALPHA);
	
if(CALF ~= 0) 
 TAN=SALF/CALF;
else
TAN=1e32;
 if (SALF < 0) 
  TAN=-TAN;
 end 
end
X=SL*TAN+VM;
end    