function decod = decodificarPorSumaProducto(recibido, H, N0vec, iteraciones)
% Se implementa el algoritmo suma-producto para realizar la decodificaci�n
% por medio del LLR (Log-likelyhood ratio)
%
%  recibdo     : Mensaje recibido
%  H           : Matriz de checkeo de paridad
%  N0vec       : Varianza del ruido
%  iteraciones : Numero de iteraciones
%
%  decod       : Mensaje decodificado 
%

[M N] = size(H);
% Se determina el LLR a priori para canal AWGN, utilizando solos los 
% valores recibidos. El signo menos es porque se mapea 0 con -1 y 1 con 1.
LLRi = (-4*recibido./N0vec)';

% Inicializaci�n de matrices
Eij = zeros(M, N);
PhiBeta = zeros(M, N);

% Se determina la matriz M inicial (mensajes de los nodos variables a los
% nodos funci�n) con el LLR a priori.
Mij = H.*repmat(LLRi, M, 1);
% Las i-�sima columna de M contiene el LLR del bit "i" en todas sus filas 
% en las que H tiene un "1"

[fil, col] = find(H);
% fil y col contienen las posiciones (dadas en fila y columna) de H en las
% que hay "1"
% Comienzan las iteraciones
for n = 1:iteraciones
   
   % Se definen las matrices alpha y beta, que contienen el signo y el 
   % m�dulo de la matriz M
   alpha = sign(Mij);   
   beta = abs(Mij);

   % Se calcula la matriz de la funci�n Phi(beta) en las posiciones
   % relevantes (en las que H tiene un "1"). Se va a aplicar la ecuaci�n
   % 2.14 del paper de Sarah Johnson
   for l = 1:length(fil)
      PhiBeta(fil(l), col(l)) = log((exp(beta(fil(l), col(l))) + 1)/...
                             (exp(beta(fil(l), col(l))) - 1));
   end
   
   % Se calcula ahora la matriz E con la f�rmula 2.14 del paper de Johnson.
   % Esta matriz simboliza los mensajes que se env�an desde los nodos
   % funci�n a los nodos variables.
   for i = 1:M
      
      % c1 tiene las columnas en las que se encuentran los elementos no
      % nulos de H para la i-�sima fila
      c1 = find(H(i, :));
      
      % Este ciclo for calcula la matriz compuesta por la sumatoria de 
      % todos los LLR de PhiBeta menos el LLR en esa componente (ver
      % ecuaci�n).
      for k = 1:length(c1)

         SumaPhiBeta = 0;
         ProductoAlpha = 1;
         
         % Sumatoria de Phi(beta)\c1(k) , es decir, excluy�ndose el valor
         % de dicha columna
         SumaPhiBeta = sum(PhiBeta(i, c1)) - PhiBeta(i, c1(k));
         
         % Esta es una validaci�n para evitar divisi�n por cero "num�rica"
         if SumaPhiBeta < 1e-20
            SumaPhiBeta = 1e-10;
         end
         
         % Se calcula Phi(sumatoria(Phi(beta)))
         PhiSumaPhiBeta = log((exp(SumaPhiBeta) + 1)/(exp(SumaPhiBeta) - 1));
      
         % Multiplication of alphaij\c1(k) (use '*' since alphaij are -1/1s)
         % Se calcula la productoria de los Alpha\c1(k) . Para evitar 
         % multiplicar por el valor de Alpha de esa componente, se vuelve a
         % multiplicar por este valor en vez de dividir por este valor, lo 
         % que puede hacerse ya que los Alpha son solamente " 1" y "-1".
         ProductoAlpha = prod(alpha(i, c1))*alpha(i, c1(k));
         
         % Se calcula entonces Eij con la ecuaci�n 2.14 del paper
         Eij(i, c1(k)) = ProductoAlpha*PhiSumaPhiBeta; %Ecuac��n 2.14
         
      end % for k
      
   end % for i

   % Se calculan ahora, bajo el mismo ciclo for, la matriz M y el LLR total
   % luego de esta iteraci�n. Con el LLR total, se puede decodificar
   % el mensaje y obtener el mensaje decodificado luego de esta iteraci�n.
   for j = 1:N

      % Find non-zero in the row
      f1 = find(H(:, j));
      
      for k = 1:length(f1)        
        
         % Se obtiene entonces la matriz Mij sumando los valores de cada
         % columna de Eij menos el valor de esa posici�n y sumando a su vez
         % el LLR a priori
         Mij(f1(k), j) = LLRi(j) + sum(Eij(f1, j)) - Eij(f1(k), j);
      
      end % for k
      
      % Se calcula el LLR de esta iteraci�n como suma de las los LLRs en
      % las columnas de Eij y el LLR a priori
      LLR = LLRi(j) + sum(Eij(f1, j));
      
      % Se decodifica con el LLR
      if LLR < 0
         decod(j) = 1;
      else
         decod(j) = 0;
      end
                       
   end % for j
   
end % for n
