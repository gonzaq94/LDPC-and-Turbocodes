function H = generarLDPC(M, N, unosPorCol)
%Funcion que crea la Matriz de Checkeo de Pararidad H para un LDPC 
%de rate 1/2(Mc....Neal)
%
%  M    : Cantidad de filas
%  N    : Cantidad de Columnas
%  unosPorCol: Cantidad de unos por columna
%
%  H    : Matriz de Checkeo de Paridad de baja densisdad                   

%Rate debe ser 1/2
unosPorFila = (N/M)*unosPorCol;

% Se genera matriz de indices de columnas desordenados en forma aleatoria
for i = 1:N
 unosEnCol(:, i) = randperm(M)';
end

% Se crea indice para lo elementos no nulos(1s) 
r = reshape(unosEnCol(1:unosPorCol, :), N*unosPorCol, 1);
tmp = repmat([1:N], unosPorCol, 1);
c = reshape(tmp, N*unosPorCol, 1);

% Se ordena en orden ascendente el indice de las filas
[r, ix] = sort(r);

% Se ordena el indice de las columnas segun el de las filas
for i = 1:N*unosPorCol
 cSort(i, :) = c(ix(i));
end

% Se crea otro indice de filas aletoriamnte distribuidos.
% Notar que ahora los indices estan totalmente desligados.
tmp = repmat([1:M], unosPorFila, 1);
r = reshape(tmp, N*unosPorCol, 1);


%Creamos H.Se elimina cualquier indice dulpicado mediante una AND
%Nota:Tener en cuenta que sparse devuleve un objeto no una matriz con lo
%cual no se realiza la AND habitual.
S = and(sparse(r, cSort, 1, M, N), ones(M, N));
H = full(S);     
      
% Se rellenan las filas que no tienen unos o bien tienen uno solo
for i = 1:M
   
   n = randperm(N);
   %Si no tiene unos se agregan dos unos en posiciones aleatorias
   if length(find(r == i)) == 0
      H(i, n(1)) = 1;
      H(i, n(2)) = 1;
   %Si tiene un solo uno se agrega otro en posicion aleatoria
   elseif length(find(r == i)) == 1
      H(i, n(1)) = 1;
   end

end 

%Eliminacion de ciclos de longitud 4
for i = 1:M
    
  %Buscacomos los cilos
  for j = (i + 1):M         
     %Miramos si hay 1s en la misma posicion en filas consecutivas
     w = and(H(i, :), H(j, :));
     c1 = find(w);
     lc = length(c1);
     if lc > 1

       
        %Al encontrar un lazo cambia 1 por 0 en la fila con mayor numero 1s
        if length(find(H(i, :))) < length(find(H(j, :)))
           %Se repiten el proceso hasta que quede solo una columna en la
           %que se repiten
           for cc = 1:lc - 1
              H(j, c1(cc)) = 0;
           end
        else
           for cc = 1:lc - 1
              H(i, c1(cc)) = 0;
           end
        end             

     end 
  end 
end 
  
end 


