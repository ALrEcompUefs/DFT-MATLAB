% Alisson Rodrigues de Carvalho Santos
% Universidade Estadual de Feira de Santana
% TEC 430 Processamento digital de sinais

clc;
clear all
close all

% Este script permite realizar a amostragem e comparação dos métodos com
% diferentes quantidades de amostras

% Gerando sinal de entrada
dt = 1/4000;
t = 0:dt:1;

% Sinal de entrada
x =  cos(2*pi*100.*t)+3*cos(2*pi*250.*t)+5*cos(2*pi*750.*t)+7*cos(2*pi*1000.*t);

% Frequência de amostragem
Fs = 2500; % em Hertz
% Período do amostragem
Ts = 1/Fs;

ts = 0:Ts:1;
% sinal amostrado
xn = cos(2*pi*100.*ts)+3*cos(2*pi*250.*ts)+5*cos(2*pi*750.*ts)+7*cos(2*pi*1000.*ts);

% Plotagem do sinal de entrada e sinal de entrada amostrado
figure('name','comparativo de sinais');
subplot(2,1,1);
plot(t,x);title('Sinal original'); grid on;
xlabel('tempo(s)'); ylabel('amplitude');
subplot(2,1,2);
stem(ts,xn);title('Sinal amostrado'); grid on;
xlabel('tempo(s)'); ylabel('amplitude');

%% Aplicação da DFT e FFT
% Numero de amostras utilizadas
N  = 32;
% Numero de amostras retiradas do sinal
% Informe um valor menor que N se quiser completar com zeros depois
N_ = N; 
xn_janelado = xn(1:N);

% Nome da figura
figura_name = 'Comparativo DFT';
titulo = [num2str(N),' amostras'];
% plot disctreto ?
discreto = false;

% Para completar a amostra com zeros mude a variavel para true e defina M
% com o total de zeros a prencher
% Considere que o complemento não dever passar de N
zero_padding = false;
% total de zeros
M =0;
if zero_padding
    xn_janelado = [xn_janelado, zeros(1,M) ];
end

% Sinal Janelado
figure;
plot(xn_janelado);grid on; title('Sinal janelado 32 amostras');
xlabel('tempo(s)'); ylabel('amplitude');

% Matrizes para armazenar resultados

% Valores da MyDFT,MyFFT e fft
Xk = zeros(3,32);
% Total de multiplicações e divisões

n_Add = zeros(2,1);
n_Mult = zeros(2,1);

% Aplicando a dft ao sinal amostrado
[Xk(1,:),n_Add(1),n_Mult(1)] = MyDFT(xn_janelado,N);
% Aplicando a fft com dizimação na frequência
[Xk(2,:),n_Add(2),n_Mult(2)] = MyFFT(xn_janelado,N);
% Aplicando função FFT nativa do matlab
Xk(3,:) = fft(xn_janelado,N);
% Função para plotagem dos espectros
plot_espectros(Xk(1,:),Xk(2,:),Xk(3,:),figura_name,titulo,N,Fs,discreto);



% Exibe no console o total de operações
fprintf('|-----MyDFT-----|-----MyFFT-----|\n');
fprintf('|--ADD  PROD----|--ADD  PROD----|\n');
fprintf('|--%d  %d----|--%d  %d----|\n',n_Add(1),n_Mult(1),n_Add(2),n_Mult(2));
fprintf('|-------------------------------|\n');



%% Função para plotagem dos gráficos
function [] = plot_espectros(xk_dft,xk_myfft,xk_fft,figura,titulo,N,Fs,discreto)
    f = (0:N-1) * Fs / N;
    figure('name',figura);
    subplot(3,1,1);
    if discreto
        stem(f,abs(xk_dft),'.');title(['MyDFT ',titulo]); grid on;
    else
        plot(f,abs(xk_dft));title(['MyDFT ',titulo]); grid on;
    end
    xlabel('frequência(Hz)');
    ylabel('Magnitude');
    subplot(3,1,2);
    if discreto
        stem(f,abs(xk_myfft),'.');title(['MyFFT ',titulo]); grid on;
    else
        plot(f,abs(xk_myfft),'.');title(['MyFFT ',titulo]); grid on;
    end
    xlabel('frequência(Hz)');
    ylabel('Magnitude');
    subplot(3,1,3);
    if discreto
        stem(f,abs(xk_fft),'.');title(['FFT ',titulo]); grid on;
    else
        plot(f,abs(xk_fft));title(['FFT ',titulo]); grid on;
    end
    xlabel('frequência(Hz)');
    ylabel('Magnitude');
end
