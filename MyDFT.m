%
%
% Alisson Rodrigues de Carvalho Santos
function [Xk,n_Add,n_Mult] = MyDFT(xn,N)
    % Vetor de coeficientes
    Xk= zeros(1,N);
    % constante 2π/N
    n_Add =0;
    n_Mult =0;
    % Primeiro for computa os k coeficientes da DFT
    for k=0:1:N-1
        % Segundo for realiza o somatorio em n do coeficiente K
        for n=0:1:N-1
            % caclula o fator wn
            wn = exp((-1i*2*pi*n*k)/N);
            Xk(k+1)= Xk(k+1)+( xn(n+1)* wn);
            n_Mult = n_Mult +1;
            n_Add = n_Add+1;
        end
        % compensa soma em n=zero
        n_Add = n_Add-1;
    end
end