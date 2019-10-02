function [c, nuevaH] = codificar(entrada, H)
% Genera un vector de paridad mediante la descomposión de H en la forma
% H=[A B], donde A=L*U. Si el mensaje codificado es u=[c m]^t, se
% debe cumplir H*u=0. Pero por la descomposición dada 
% H*u=[A B]*[c m]^t=A*c^t+B*m^t=0, lo que implica que A*c^t=B*m^t.
% Definiendo z=B*m^t, A*c^t=L*U**c^t=z. Como L y U son matrices LU, el
% despeje puede hacerse de manera fácil.
%
%  entrada : Mensaje binario a codificar (0/1)
%  H       : Matriz de checkeo de paridad
%           
%  c       : Bits de paridad
%  nuevaH    : H con las columnas alteradas (debe usarse para la
%            decodificación en lugar de la H original, pues las ecuaciones de
%            checkeo de paridad cambiaron).
%

% Inicializaciones
[M, N] = size(H);

F = H;

L = zeros(M, N - M);
U = zeros(M, N - M);

% Se procede a hacer permutaciones de filas y columnas que pongan 1s en los
% elementos de la diagonal de la submatriz A (las primeras N-M columnas)
for i = 1:M
    
         % Encuentra todas las posiciones donde hay elementos no nulos. Las
         % columnas anteriores a la i-ésima no las toco para no arruinar la
         % descomposición ya hecha.
         [fils, cols] = find(F(:, i:end));
         % Se calculan los "pesos" de las filas y columnas, esto es, se
         % calcula la cantidad de 1s que hay en cada columna y en la fila en
         % la que se quiere poner en 1 en la diagonal
         PesoCols = sum(F(:, i:end), 1) - 1;
         PesoFils = sum(F(i, :), 2) - 1;
         
         % Se buscan las posiciones o indices en el vector fils en el que
         % estan los 1s que se encuentran sobre la fila en la que se quiere
         % poner un 1 en la diagonal. Esto se hace para obtener en que
         % columnas estan esos 1s y poder hacer la permutación.
         IndiceFilas = find(fils == i);
            
         % Se encuentra el mínimo producto entre los pesos de las columnas
         % y los pesos de las filas. De todas las columnas que tienen un 1
         % en la fila deseada, se elegirá la que minimice el producto de
         % los pesos, para así tener una menor cantidad de 1s fuera de la
         % diagonal.
         [minimo, IndiceMin] = min(PesoCols(cols(IndiceFilas))*PesoFils);
         % Se obtiene así la columna elegida para realizar la permutación.
         % Como solo se "miraban" las columnas de indice mayor a "i", se
         % debe sumar i-1 al número de columna.
         ColElegida = cols(IndiceFilas(IndiceMin)) + (i - 1);

   % Permutación de columnas
   temp1 = F(:, i);
   temp2 = H(:, i);
   F(:, i) = F(:, ColElegida);
   H(:, i) = H(:, ColElegida);
   F(:, ColElegida) = temp1;
   H(:, ColElegida) = temp2;
                     
   % Lleno las matrices L y U
   L(i:end, i) = F(i:end, i);
   U(1:i, i) = F(1:i, i);
         
   % En la última fila no es necesario realizar ninguna permutación.
   if i < M           
            
      % Encuentra las filas que, en la columna i, tiene 1
      [f2, c2] = find(F((i + 1):end, i));          
      % Se le suma la fila actual (que se sabe que tiene un 1 en la columna
      % i) a las filas encontradas para eliminar esos 1s. De esto modo
      % aseguro que la submatriz de F de las primeras M-N=N/2 columnas
      % tiene todos 0s por debajo de la diagonal principal.
      F((i + f2), :) = mod(F((i + f2), :) + repmat(F(i, :), length(f2), 1), 2);
                                                           
   end % if
         
end % for i

% Se puede realizar entonces la descomposión LU de A (que es inversible 
% pues no tiene ceros en la diagonal).
A=L*U;
B=H(:, (N - M) + 1:end);

% Conlo que puede obtenerse el vector checkeo de paridad
z = mod(B*entrada, 2);
c = mod(U\(L\z), 2); 

% La matriz H que se utiliza para la decodificación debe modificarse
% al haberse alterado el orden de las columnas)
nuevaH = H;
