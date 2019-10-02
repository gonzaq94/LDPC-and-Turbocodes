function [ber1_tc,ber2_tc]=berTurbocodigo
%Funcion que simula la transimicion de mensajes codificados mediante un
%Turbo Codigo de rate 1/2 con modulacion BPSK sobre un canal Gilbert-Elliot y AWGN.
%Se grafica la relacion bit error rate en funcion de la SNR para los
%diferentes canales.

clc
clear all

%Cantidad de fit utilies por mensaje
M=500;

%Diferencia entre las varianzas del estado "bueno" y "malo" en dB
deltaN=3; 

%Probabilidad de transicion del estado "malo" al "bueno"
Pgb=0.3;

%Probabilidad de transicion del estado "bueno" al "malo"
Pbg=0.2;

frame=100;

% EbN0 en dB
EbN0 = [0 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3];

ber1=zeros(size(EbN0));
ber2=zeros(size(EbN0));
ber3=zeros(size(EbN0));
ber4=zeros(size(EbN0));

interleaver = randperm(M);

TurboCodificador = comm.TurboEncoder('TrellisStructure',poly2trellis(4, ...
    [13 15],13),'InterleaverIndices',interleaver);

TurboDecodificador = comm.TurboDecoder('TrellisStructure',poly2trellis(4, ...
    [13 15],13),'InterleaverIndices',interleaver, ...
    'NumIterations',1);

for i=1:length(EbN0)
    
    mensaje = round(rand(M, frame));
    
    for j=1:frame
       
       %Codificacion
       msjCodificado = step(TurboCodificador,mensaje(:,j));
       
       %Modulacion
       msjModulado = 2*msjCodificado - 1;
       
       %Simulacion Canales
       %(1)AWGN
       N0 = 1/(exp(EbN0(i)*log(10)/10));
       msjTransmitido = msjModulado + sqrt(N0)*randn(size(msjModulado));
       
       %(2) Gilbert Elliot
        Nb=N0;
        Ng=1/(exp((EbN0(i)+deltaN)*log(10)/10));
        varianza=gilbertElliot(length(msjModulado),Nb,Ng,Pbg,Pgb);
        NormalEstandar=randn(size(msjModulado));
        ruido=sqrt(varianza).*NormalEstandar;
        msjTransmitidoGE=msjModulado + ruido;
        
       %Decodificacion
       msjDecodificado1 = step(TurboDecodificador,msjTransmitido);
        
       msjDecodificado2 = step(TurboDecodificador,msjTransmitidoGE);

            
       %Calculo de Bit Error Rate     
       [num1, rat1] = biterr(msjDecodificado1, mensaje(:,j));
       ber1(i) = (ber1(i) + rat1);
       
       [num2, rat2] = biterr(msjDecodificado2, mensaje(:,j));
       ber2(i) = (ber2(i) + rat2);

    end
    
    ber1(i)=ber1(i)/frame;
    ber2(i)=ber2(i)/frame;

end

ber1_tc=ber1;
ber2_tc=ber2;
       

