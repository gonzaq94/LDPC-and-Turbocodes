function [ber1_ldpc,ber2_ldpc,ber3_ldpc,ber4_ldpc]=berLDPC
%Funcion que simula la transimicion de mensajes codificados mediante un LDPC
%de rate 1/2 con modulacion BPSK sobre un canal Gilbert-Elliot y AWGN.
%Se devuelve la relacion bit error rate en funcion de la SNR para diferentes
%decodificaciones.

clc; 
clear all;

% Cantidad de bits del msj codificado
N = 1000;

% Cantidad de bits del msj util
M=N/2;

%Diferencia entre las varianzas del estado "bueno" y "malo" en dB
deltaN=3; 

%Probabilidad de transicion del estado "malo" al "bueno"
Pgb=0.3;

%Probabilidad de transicion del estado "bueno" al "malo"
Pbg=0.2;

% Numero de 1's por columna
unosPorCol = 3;

% EbN0 en dB
EbN0 = [0 0.5 1 1.5 2 2.5 3];

% Número de iteraciones;
k = 5;

%Cantidad de mensajes a transmitir
frame = 10;

% Se crea la matriz de Checkeo de Paridad de LDPC
H = generarLDPC(M, N, unosPorCol);

for i = 1:length(EbN0)
   
   ber1(i) = 0;
   ber2(i) = 0;
   ber3(i) = 0;
   ber4(i) = 0;
   
   % Se genera el frame de msjs aleatorios de 0s y 1s a transmitir
   mensaje = round(rand(M, frame));
   
   for j = 1:frame
      
      % Codificación del mensaje
      [paridad, nuevaH] = codificar(mensaje(:, j), H);
      msjCodificado = [paridad; mensaje(:, j)];

      % Modulación de la secuencia binaria para transmitir una señal
      % analógica
      msjModulado = 2*msjCodificado - 1;
      % Canal AWGN
      N0 = 1/(exp(EbN0(i)*log(10)/10));
      msjTransmitido = msjModulado + sqrt(N0)*randn(size(msjModulado));
      
      % Canal Gillbert Elliot
      Nb=N0;
      Ng=1/(exp((EbN0(i)+deltaN)*log(10)/10));
      varianza=gilbertElliot(N,Nb,Ng,Pbg,Pgb);
      NormalEstandar=randn(N,1);
      ruido=sqrt(varianza).*NormalEstandar;
      msjTransmitidoGE=msjModulado + ruido;

      % Decodificación
      
      N0Vec=N0*ones(N,1);
      msjDecodificado1 = decodificarPorSumaProducto(msjTransmitido, nuevaH, N0Vec, k);
      
      NbVec=Nb*ones(N,1);
      msjDecodificado2 = decodificarPorSumaProducto(msjTransmitidoGE, nuevaH,NbVec, k);
      
      NgVec=Ng*ones(N,1);
      msjDecodificado3 = decodificarPorSumaProducto(msjTransmitidoGE, nuevaH,NgVec, k);
      
      VectorVarianzas=gilbertElliot(N,Nb,Ng,Pbg,Pgb);
      msjDecodificado4 = decodificarPorSumaProducto(msjTransmitidoGE, nuevaH,VectorVarianzas, k);
     
   
      % Cálculo del BER (se incluyen los bits de paridad)
      [num1, rat1] = biterr(msjDecodificado1', msjCodificado);
      ber1(i) = (ber1(i) + rat1);
      
      [num2, rat2] = biterr(msjDecodificado2', msjCodificado);
      ber2(i) = (ber2(i) + rat2);
      
      [num3, rat3] = biterr(msjDecodificado3', msjCodificado);
      ber3(i) = (ber3(i) + rat3);
      
      [num4, rat4] = biterr(msjDecodificado4', msjCodificado);
      ber4(i) = (ber4(i) + rat4);
      
   end % for j
   
   % Se toman los promedios de los BER
   ber1(i) = ber1(i)/frame;
   ber2(i) = ber2(i)/frame;
   ber3(i) = ber3(i)/frame;
   ber4(i) = ber4(i)/frame;
   
end % for i

ber1_ldpc=ber1;
ber2_ldpc=ber2;
ber3_ldpc=ber3;
ber4_ldpc=ber4;