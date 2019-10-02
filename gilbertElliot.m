function [varianza] = gilbertElliot (M,Nb,Ng,Pbg,Pgb)

% Se implementa el canal con memoria Gillbert Elliot. El ruido en este tipo
% de canal es una cadena de Markov con dos estados, uno bueno y uno malo,
% en el que las varianzas del ruido gaussiano son distintas.
%
%  M     : Largo del vector de ruido
%  Nb    : Varianza en el estado "malo"
%  Ng    : Varianza en el estado "bueno"
%  Pbg     : Probabilidad de transición del estado "bueno" al "malo"
%  Pgb     : Probabilidad de transición del estado "malo" al "bueno"
%
%  varianza : Vector de varianzas 


varianza=zeros(M,1);

varianza(1)=Ng;

for i=2:length(varianza)
    Naleatorio=rand;
    if varianza(i-1)==Ng
        if Naleatorio>Pbg 
            varianza(i)=Ng;
        else
            varianza(i)=Nb;
        end
    else 
        if Naleatorio>Pgb 
            varianza(i)=Nb;
        else
            varianza(i)=Ng;
        end
    end
end

