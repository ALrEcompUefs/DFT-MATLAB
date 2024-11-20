% Alisson Rodrigues de Carvalho Santos
% TEC 430 Processamento digital de sinais

% Segunda Avaliação
% Implementação da DFT e FFT com dizimação na frequência
clc;
clear all;
close all;

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

figure('name','comparativo de sinais');
subplot(2,1,1);
plot(t,x);title('Sinal original'); grid on;
xlabel('tempo(s)'); ylabel('amplitude');
subplot(2,1,2);
stem(ts,xn);title('Sinal amostrado'); grid on;
xlabel('tempo(s)'); ylabel('amplitude');

% Para cada requisito da atividade avaliativa foi feita uma seção neste
% codigo que aplica as funções MyDFT,MyFFT e a fft do matlab
% Foram implementadas as função plot_espectro e plot_espectro que realizam
% a plotagem das figuras dos espectros resultantes das funções

%%  Janelamento de 32 amostras ponto a
%
% Aplicando o janelamento no sinal para obter as amostras
N= 32;
xn_janelado = xn(1:32);

figure;
plot(xn_janelado);grid on; title('Sinal janelado 32 amostras');
xlabel('tempo(s)'); ylabel('amplitude');

% matriz para armazenar resultados
Xk = zeros(3,32);
n_Add = zeros(2,1);
n_Mult = zeros(2,1);
% Aplicando a dft ao sinal amostrado
[Xk(1,:),n_Add(1),n_Mult(1)] = MyDFT(xn_janelado,N);
% Aplicando a fft com dizimação na frequência
[Xk(2,:),n_Add(2),n_Mult(2)] = MyFFT(xn_janelado,N);
% Aplicando função FFT nativa do matlab
Xk(3,:) = fft(xn_janelado,N);
% Função para plotagem dos espectros
plot_espectros(Xk(1,:),Xk(2,:),Xk(3,:),'FFT 32 amostras',N,Fs);

%% Completando 64 amostras com zeros a direita ponto b
N= 64;
xn_zero_padding = [xn_janelado, zeros(1,32) ];

figure;
plot(xn_zero_padding);grid on; title('Zero padding 64 amostras');
xlabel('tempo(s)'); ylabel('amplitude');xlim([0,N]);

% matriz para armazenar resultados
Xk = zeros(3,64);
n_Add = zeros(2,1);
n_Mult = zeros(2,1);
% Aplicando a dft ao sinal amostrado
[Xk(1,:),n_Add(1),n_Mult(1)] = MyDFT(xn_zero_padding,N);
% Aplicando a fft com dizimação na frequência
[Xk(2,:),n_Add(2),n_Mult(2)] = MyFFT(xn_zero_padding,N);
% Aplicando função FFT nativa do matlab
Xk(3,:) = fft(xn_zero_padding,N);
% Função para plotagem dos espectros
plot_espectros2(Xk(1,:),Xk(2,:),Xk(3,:),'FFT 64 amostras',N,Fs);

%% Janelamento 64 amostras ponto c
%
% Aplicando o janelamento no sinal para obter as amostras
N= 64;
xn_janelado = xn(1:N);

figure;
plot(xn_janelado);grid on; title('Sinal janelado 64 amostras');
xlabel('tempo(s)'); ylabel('amplitude');xlim([0,N])

% matriz para armazenar resultados
Xk = zeros(3,N);
n_Add = zeros(2,1);
n_Mult = zeros(2,1);
% Aplicando a dft ao sinal amostrado
[Xk(1,:),n_Add(1),n_Mult(1)] = MyDFT(xn_janelado,N);
% Aplicando a fft com dizimação na frequência
[Xk(2,:),n_Add(2),n_Mult(2)] = MyFFT(xn_janelado,N);
% Aplicando função FFT nativa do matlab
Xk(3,:) = fft(xn_janelado,N);
% Função para plotagem dos espectros
plot_espectros2(Xk(1,:),Xk(2,:),Xk(3,:),'FFT 64 amostras',N,Fs);

%% Completando 128 amostras com zeros a direita ponto d
N= 128;
xn_zero_padding = [xn_janelado, zeros(1,64) ];

figure;
plot(xn_zero_padding);grid on; title('Zero padding 128 amostras');
xlabel('tempo(s)'); ylabel('amplitude');xlim([0,N]);

% matriz para armazenar resultados
Xk = zeros(3,N);
n_Add = zeros(2,1);
n_Mult = zeros(2,1);
% Aplicando a dft ao sinal amostrado
[Xk(1,:),n_Add(1),n_Mult(1)] = MyDFT(xn_zero_padding,N);
% Aplicando a fft com dizimação na frequência
[Xk(2,:),n_Add(2),n_Mult(2)] = MyFFT(xn_zero_padding,N);
% Aplicando função FFT nativa do matlab
Xk(3,:) = fft(xn_zero_padding,N);
% Função para plotagem dos espectros
plot_espectros2(Xk(1,:),Xk(2,:),Xk(3,:),'FFT 128 amostras',N,Fs);

%% Janelamento 128 amostras
%
% Aplicando o janelamento no sinal para obter as amostras
N= 128;
xn_janelado = xn(1:N);

figure;
plot(xn_janelado);grid on; title('Sinal janelado 128 amostras');
xlabel('tempo(s)'); ylabel('amplitude');xlim([0,N])

% matriz para armazenar resultados
Xk = zeros(3,N);
n_Add = zeros(2,1);
n_Mult = zeros(2,1);
% Aplicando a dft ao sinal amostrado
[Xk(1,:),n_Add(1),n_Mult(1)] = MyDFT(xn_janelado,N);
% Aplicando a fft com dizimação na frequência
[Xk(2,:),n_Add(2),n_Mult(2)] = MyFFT(xn_janelado,N);
% Aplicando função FFT nativa do matlab
Xk(3,:) = fft(xn_janelado,N);
% Função para plotagem dos espectros
plot_espectros2(Xk(1,:),Xk(2,:),Xk(3,:),'FFT 128 amostras',N,Fs);

%% Completando 256 amostras com zeros a direita
N= 256;
xn_zero_padding = [xn_janelado, zeros(1,128) ];

figure;
plot(xn_zero_padding);grid on; title('Zero padding 256 amostras');
xlabel('tempo(s)'); ylabel('amplitude');xlim([0,N]);

% matriz para armazenar resultados
Xk = zeros(3,N);
n_Add = zeros(2,1);
n_Mult = zeros(2,1);
% Aplicando a dft ao sinal amostrado
[Xk(1,:),n_Add(1),n_Mult(1)] = MyDFT(xn_zero_padding,N);
% Aplicando a fft com dizimação na frequência
[Xk(2,:),n_Add(2),n_Mult(2)] = MyFFT(xn_zero_padding,N);
% Aplicando função FFT nativa do matlab
Xk(3,:) = fft(xn_zero_padding,N);
% Função para plotagem dos espectros
plot_espectros2(Xk(1,:),Xk(2,:),Xk(3,:),'FFT 256 amostras',N,Fs);

%% Janelamento 256 amostras
%
% Aplicando o janelamento no sinal para obter as amostras
N= 256;
xn_janelado = xn(1:N);

figure;
plot(xn_janelado);grid on; title('Sinal janelado 128 amostras');
xlabel('tempo(s)'); ylabel('amplitude');xlim([0,N])

% matriz para armazenar resultados
Xk = zeros(3,N);
n_Add = zeros(2,1);
n_Mult = zeros(2,1);
% Aplicando a dft ao sinal amostrado
[Xk(1,:),n_Add(1),n_Mult(1)] = MyDFT(xn_janelado,N);
% Aplicando a fft com dizimação na frequência
[Xk(2,:),n_Add(2),n_Mult(2)] = MyFFT(xn_janelado,N);
% Aplicando função FFT nativa do matlab
Xk(3,:) = fft(xn_janelado,N);
% Função para plotagem dos espectros
plot_espectros2(Xk(1,:),Xk(2,:),Xk(3,:),'FFT 128 amostras',N,Fs);

%% Completando 512 amostras com zeros a direita
N= 512;
xn_zero_padding = [xn_janelado, zeros(1,256) ];

figure;
plot(xn_zero_padding);grid on; title('Zero padding 512 amostras');
xlabel('tempo(s)'); ylabel('amplitude');xlim([0,N]);

% matriz para armazenar resultados
Xk = zeros(3,N);
n_Add = zeros(2,1);
n_Mult = zeros(2,1);
% Aplicando a dft ao sinal amostrado
[Xk(1,:),n_Add(1),n_Mult(1)] = MyDFT(xn_zero_padding,N);
% Aplicando a fft com dizimação na frequência
[Xk(2,:),n_Add(2),n_Mult(2)] = MyFFT(xn_zero_padding,N);
% Aplicando função FFT nativa do matlab
Xk(3,:) = fft(xn_zero_padding,N);
% Função para plotagem dos espectros
plot_espectros2(Xk(1,:),Xk(2,:),Xk(3,:),'FFT 512 amostras',N,Fs);

%% Janelamento 512 amostras
%
% Aplicando o janelamento no sinal para obter as amostras
N= 512;
xn_janelado = xn(1:N);

figure;
plot(xn_janelado);grid on; title('Sinal janelado 512 amostras');
xlabel('tempo(s)'); ylabel('amplitude');xlim([0,N])

% matriz para armazenar resultados
Xk = zeros(3,N);
n_Add = zeros(2,1);
n_Mult = zeros(2,1);
% Aplicando a dft ao sinal amostrado
[Xk(1,:),n_Add(1),n_Mult(1)] = MyDFT(xn_janelado,N);
% Aplicando a fft com dizimação na frequência
[Xk(2,:),n_Add(2),n_Mult(2)] = MyFFT(xn_janelado,N);
% Aplicando função FFT nativa do matlab
Xk(3,:) = fft(xn_janelado,N);
% Função para plotagem dos espectros
plot_espectros2(Xk(1,:),Xk(2,:),Xk(3,:),'FFT 512 amostras',N,Fs);

%% Completando 1024 amostras com zeros a direita
N= 1024;
xn_zero_padding = [xn_janelado, zeros(1,512) ];

figure;
plot(xn_zero_padding);grid on; title('Zero padding 1024 amostras');
xlabel('tempo(s)'); ylabel('amplitude');xlim([0,N]);

% matriz para armazenar resultados
Xk = zeros(3,N);
n_Add = zeros(2,1);
n_Mult = zeros(2,1);
% Aplicando a dft ao sinal amostrado
[Xk(1,:),n_Add(1),n_Mult(1)] = MyDFT(xn_zero_padding,N);
% Aplicando a fft com dizimação na frequência
[Xk(2,:),n_Add(2),n_Mult(2)] = MyFFT(xn_zero_padding,N);
% Aplicando função FFT nativa do matlab
Xk(3,:) = fft(xn_zero_padding,N);
% Função para plotagem dos espectros
plot_espectros(Xk(1,:),Xk(2,:),Xk(3,:),'FFT 1024 amostras',N,Fs);

%% Janelamento 1024 amostras
%
% Aplicando o janelamento no sinal para obter as amostras
N= 1024;
xn_janelado = xn(1:N);

figure;
plot(xn_janelado);grid on; title('Sinal janelado 1024 amostras');
xlabel('tempo(s)'); ylabel('amplitude');xlim([0,N])

% matriz para armazenar resultados
Xk = zeros(3,N);
n_Add = zeros(2,1);
n_Mult = zeros(2,1);
% Aplicando a dft ao sinal amostrado
[Xk(1,:),n_Add(1),n_Mult(1)] = MyDFT(xn_janelado,N);
% Aplicando a fft com dizimação na frequência
[Xk(2,:),n_Add(2),n_Mult(2)] = MyFFT(xn_janelado,N);
% Aplicando função FFT nativa do matlab
Xk(3,:) = fft(xn_janelado,N);
% Função para plotagem dos espectros
plot_espectros2(Xk(1,:),Xk(2,:),Xk(3,:),'FFT 1024 amostras',N,Fs);

%% Completando 2048 amostras com zeros a direita
N= 2048;
xn_zero_padding = [xn_janelado, zeros(1,1024) ];

figure;
plot(xn_zero_padding);grid on; title('Zero padding 2048 amostras');
xlabel('tempo(s)'); ylabel('amplitude');xlim([0,N]);

% matriz para armazenar resultados
Xk = zeros(3,N);
n_Add = zeros(2,1);
n_Mult = zeros(2,1);
% Aplicando a dft ao sinal amostrado
[Xk(1,:),n_Add(1),n_Mult(1)] = MyDFT(xn_zero_padding,N);
% Aplicando a fft com dizimação na frequência
[Xk(2,:),n_Add(2),n_Mult(2)] = MyFFT(xn_zero_padding,N);
% Aplicando função FFT nativa do matlab
Xk(3,:) = fft(xn_zero_padding,N);
% Função para plotagem dos espectros
plot_espectros2(Xk(1,:),Xk(2,:),Xk(3,:),'FFT 2048 amostras',N,Fs);

%%
function [] = plot_espectros(xk_dft,xk_myfft,xk_fft,titulo,N,Fs)
    f = (0:N-1) * Fs / N;
    figure('name',titulo);
    subplot(3,1,1);
    plot(f,abs(xk_dft));title('MyDFT 2048 amostras'); grid on;
    xlabel('frequência(Hz)');
    ylabel('Magnitude');
    subplot(3,1,2);
    plot(f,abs(xk_myfft));title('MyFFT 2048 amostras'); grid on;
    xlabel('frequência(Hz)');
    ylabel('Magnitude');
    subplot(3,1,3);
    plot(f,abs(xk_fft));title('FFT 2048 amostras'); grid on;
    xlabel('frequência(Hz)');
    ylabel('Magnitude');
end

function [] = plot_espectros2(xk_dft,xk_myfft,xk_fft,titulo,N,Fs)
    f = (0:N-1) * Fs / N;
    figure('name',titulo);
    subplot(3,1,1);
    stem(f,abs(xk_dft));title('MyDFT 2048 amostras'); grid on;
    xlabel('frequência(Hz)');
    ylabel('Magnitude');
    subplot(3,1,2);
    stem(f,abs(xk_myfft));title('MyFFT 2048 amostras'); grid on;
    xlabel('frequência(Hz)');
    ylabel('Magnitude');
    subplot(3,1,3);
    stem(f,abs(xk_fft));title('FFT 2048 amostras'); grid on;
    xlabel('frequência(Hz)');
    ylabel('Magnitude');
end