%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Trabajo Practico Final 
%Asignatura:Teoria de la Informacion y Decodificacion
%Autores:Gonzalo Quintana      Pad.:96689
%        Nicolás Cassia                 Pad.:96031
%Fecha: 1ro cuat 2017
%Tema: Comparacion del uso de LDPC y Turbo Codigos en la codificacion sobre
%       un canal GilbertElliot.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EbN0_ldpc= [0 0.5 1 1.5 2 2.5  3];
EbN0_tc = [0 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3];
%Bit error rate usando LDPC
[ber1_ldpc,ber2_ldpc,ber3_ldpc,ber4_ldpc]=berLDPC;

%Bit error rate usando Turbocodigos
[ber1_tc,ber2_tc]=berTurbocodigo;

%Grafico para analisis
% Gráfico del BER para LDPC en canales AWGN y Gilbert-Elliot
figure
semilogy(EbN0_ldpc, ber1_ldpc, 'o-');
hold on;
semilogy(EbN0_ldpc, ber2_ldpc, 'o--');
hold on;
semilogy(EbN0_ldpc, ber3_ldpc, 'o--');
hold on;
semilogy(EbN0_ldpc, ber4_ldpc, 'o--');
grid on;
hold off;
xlabel('SNR (dB)');
ylabel('BER');
legend('AWGN','GE decoficado con Nb','GE decoficado con Ng','GE decodificado estimando el canal');

print('BER LDPC.png','-dpng');

% Gráfico del BER para Turbo codigo en canales AWGN y Gilbert-Elliot
figure
semilogy(EbN0_tc,ber1_tc,'o-');
hold on
semilogy(EbN0_tc, ber2_tc, 'o--');
grid on;
xlabel('SNR (dB)');
ylabel('BER');
legend('AWGN','Gilbert-Elliot');

print('BER Turbocode.png','-dpng');

% Gráfico BER en AWGN para Turbo y LDPC
figure
semilogy(EbN0_tc,ber1_tc,'o-');
hold on
semilogy(EbN0_ldpc, ber1_ldpc, 'o--');
grid on;
xlabel('SNR (dB)');
ylabel('BER');
legend('AWGN Turbo','AWGN LDPC');

print('AWGN.png','-dpng');

% Gráfico BER en Gilbert-Elliot para Turbo yLDPC 
figure
semilogy(EbN0_tc, ber2_tc, 'o-');
hold on;
semilogy(EbN0_ldpc, ber2_ldpc, 'o--');
hold on;
semilogy(EbN0_ldpc, ber3_ldpc, 'o--');
hold on;
semilogy(EbN0_ldpc, ber4_ldpc, 'o--');
grid on;
hold off;
xlabel('SNR (dB)');
ylabel('BER');
legend('GE Turbo','GE decoficado con Nb','GE decoficado con Ng','GE decodificado estimando el canal');

print('GE.png','-dpng');
