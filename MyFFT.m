function [Xk, n_Add, n_Mult] = MyFFT(xn, N)
    % A função é implementada recursivamente
    % A condição de parada é um único elemento no vetor
    if N == 1
        % Retorna o valor da amostra atual
        Xk = xn;
        % Sem operações
        n_Add = 0;
        n_Mult = 0;
    else
        % Divide os dados em pares e ímpares
        xp = xn(1:2:end); % Elementos pares
        xi = xn(2:2:end); % Elementos ímpares
        % Chamada recursiva
        [Xpk, n_Add_Par, n_Mult_Par] = MyFFT(xp, length(xp));
        [Xik, n_Add_Impar, n_Mult_Impar] = MyFFT(xi, length(xi));
        % Calcula o fator de peso
        w = exp(-2i * pi * (0:(N/2 - 1)) / N);
        % Combina as FFTs
        Xk = [Xpk + w.*Xik, Xpk - w.*Xik];
        % Atualiza os contadores
        n_Add = n_Add_Par + n_Add_Impar + N;
        n_Mult = n_Mult_Par + n_Mult_Impar + N / 2;
    end
end
