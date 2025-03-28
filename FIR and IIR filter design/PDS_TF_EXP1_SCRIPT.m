%% Setup inicial
clear variables;
close all;
clc;

Tc = 5; % Opcion de elegir entre 4 a 6 segundos
fs = 8000;
r = 24;
ID = 1;
%% Captura de audio y generacion de señal combinada

% Captura de frase:"Procesamiento Digital de señales 2024"
disp("Grabando x1")
audioobj1 = grabar_audio(Tc, fs, r, ID);
x1 = getaudiodata(audioobj1)';

% Grabar la señal x1 en archivo de audio
%audiowrite('Audios\x1.wav', x1, fs);

% Captura de silbido tono constante
disp("Grabando x2")
audioobj2 = grabar_audio(Tc, fs, r, ID);
x2 = getaudiodata(audioobj2)';

% Grabar la señal x2 en archivo de audio
%audiowrite('Audios\x2.wav', x2, fs);

% Crear señal combinada
x = x1 + x2;

% Grabar la señal x en archivo de audio
%audiowrite('Audios\x.wav', x, fs);

%% Graficos magnitude spectrum y spectrogram

% Graficar Magnitude spectrum
Nx = length(x); % Numero de muestras de la señal
fx = (0:Nx/2)' * fs/Nx; % Dominio de frecuencias

% Transformada de fourier
x1_fft = fft(x1);
x2_fft = fft(x2);
x_fft = fft(x);

% Frecuencias positivas
x1_fft = x1_fft(1:Nx/2+1);
x2_fft = x2_fft(1:Nx/2+1);
x_fft = x_fft(1:Nx/2+1);

% Grafico Magnitude spectrum
figure(1)

% Magnitude spectrum x1
subplot(3, 2, 1)
plot(fx*1e-3, 20*log10(abs(x1_fft)), 'b', 'linewidth', 2)
xlabel('$f$ (kHz)', 'interpreter','latex')
ylabel('$20\log_{10}|X(f)|$ (dB)', 'interpreter','latex')
title('Audio Frase');
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-2:3), 'xticklabel', 2.^(-2:3))
xlim([.1 4])
grid on
box on

% Magnitude spectrum x2
subplot(3, 2, 3)
plot(fx*1e-3, 20*log10(abs(x2_fft)), 'b', 'linewidth', 2)
xlabel('$f$ (kHz)', 'interpreter','latex')
ylabel('$20\log_{10}|X(f)|$ (dB)', 'interpreter','latex')
title('Audio Silbido');
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-2:4), 'xticklabel', 2.^(-2:4))
xlim([.1 4])
grid on
box on

% Magnitude spectrum x
subplot(3, 2, 5)
plot(fx*1e-3, 20*log10(abs(x_fft)), 'b', 'linewidth', 2)
xlabel('$f$ (kHz)', 'interpreter','latex')
ylabel('$20\log_{10}|X(f)|$ (dB)', 'interpreter','latex')
title('Audio Mezclado');
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-2:4), 'xticklabel', 2.^(-2:4))
xlim([.1 4])
grid on
box on

% Grafico espectrograma 

%Espectrograma x1
subplot(3, 2, 2)
spectrogram(x1, [], [], [], fs, 'yaxis')
%xlabel("Tiempo (s)")

% Espectrograma x2
subplot(3, 2, 4)
spectrogram(x2, [], [], [], fs, 'yaxis')

% Espectrograma x
subplot(3, 2, 6)
spectrogram(x, [], [], [], fs, 'yaxis')

%% Diseño de filtro: Rechaza banda

% Indicar de la frecuencia del tono f_inf, f_sup y bw
f_inf = 1280; % En KHz
f_sup = 1450; % En KHz
bw = f_sup - f_inf ;

% Conversion a frecuencia normalizada
theta_inf = (f_inf/fs)*2*pi;
theta_sup = (f_sup/fs)*2*pi;
theta_bw = theta_sup - theta_inf;

% Establezco el orden del filtro % CORREGIR PARA CASO ORDEN PAR
N = 81;
n = 0:N;

% Calculo filtro prototipo
ha = (sin(theta_inf*(n-N/2)))./(pi.*(n-N/2)); % Pasa bajas
hb = ((-1).^(n)).*(sin((pi-theta_sup).*(n-N/2)))./(pi.*(n-N/2)); % Pasa altas
w = blackman(N+1)';
hp = ha + hb;
h = hp .* w;

% Aplico el filtro a x
y_FIR = filter(h, 1, x);
%audiowrite("Audios\x_FIR_filtered.wav", y_FIR, fs)

% Datos en frecuencia
[HP, FREC] = freqz(hp, 1, 8192, 2*pi); % Filtro prototipo sin aventanear
[W, ~] = freqz(w, 1, 8192, 2*pi); % Ventana
[H, Hf] = freqz(h, 1, 8192, 2*pi); % Filtro aventaneado

%% Graficos de filtro FIR

% Graficos filtro en tiempo
figure(2)

subplot(4, 1, 1)
plot(h); % Filtro aventaneado
title("Respuesta impulsiva h(n)")
grid on

% Graficos filtro en frecuencia
subplot(4, 1, 2)
plot(FREC, abs(H));
title("Respuesta en frecuencia H(n): Tanto por uno")
grid on

% Magnitude spectrum H
subplot(4, 1, 3)
plot(Hf, 20*log10(abs(H)));
title("H(e^(j*theta))");
grid on

% Fase de respuesta en frecuencia
subplot(4, 1, 4)
plot(Hf, angle(H));
title("H(e^(j*theta))");
grid on

% Filtro prototipo y ventana en el tiempo
%subplot(3, 2, 1)
%plot(hp); % Filtro prototipo
%grid on
%
%subplot(3, 2, 3)
%plot(w); % Ventana
%grid on

% Graficos tanto por uno prototipo y ventana en frecuencia
% subplot(3, 2, 2)
% plot(FREC, abs(HP)); % Filtro prototipo
% title("HP");
% grid on
% 
% subplot(3, 2, 4)
% plot(FREC, abs(W)); % Ventana
% title("W");
% grid on



%% Filtro IIR

p = butterworth_coeff_IIR(13, tan(0.9425/2));
hz = butterworth_syms_IIR(p, tan(0.9425/2), tan(0.9425/2), tan(1.2566/2));

% Otros filtros
% Orden 8, theta = 0.9425
%B = [3.5849e-04 0.0029 0.01 0.0201 0.0251 0.0201 0.01 0.0029 3.5849e-04];
%A = [1 -3.1845 5.1826 -5.2158 3.4952 -1.5728 0.4607 -0.0798 0.0062];
%
%Orden 8, theta_inf = pi/4
%B = [1.0791e-04 8.6329e-04 0.0030 0.0060 0.0076 0.0060 0.0030 8.6329e-04 1.0791e-04];
%A = [1 - 3.9838 7.5362 -8.5998 6.4002 -3.1560 1.0017 -0.1863 0.0155];
%
% Orden 8, theta_inf = pi/4, wp1 = 1200*2*pi, wp2 = 1600*2*pi
%B = [1 -8 28 -56 70 -56 28 -8 1];
%A = [1 -8 28 -56 70 -56 28 -8 1];
%
% Orden 2, theta = pi/4
%B = [0.0976 0.1953 0.0976];
%A = [1 -0.9428 0.3333];
%
% Extraido manualmente de hz
% Orden 2, tp = pi/4, wp1 = 0.9425, wp2 = 1.2566
%B = [1.2482 -2.2950 3.5514 -2.2950 1.2482 ];
%A = [1 -2.0380 3.4888 -2.5520 1.5591];
%
% N3 tc = pi/4, wl=1200*2*pi,  wh=1600*2*pi
%B = [1 3 -3 1];
%A = [1 3 -3 1];
% Pasabajos tc = pi/4 wp1 = 0.9425, wp2 = 1.2566
%
%B = [0.0317 0.0951 0.0951 0.0317];
%A = [-1 1.459 -0.9104 0.1978];
% Pasabanda tc = tan(0.9425/2), N=8, wpl = tan(0.9425/2), wph=tan(1.2566/2)
%B = [1.7463 -22.9073 145.4332 -591.4661 1.7213e+03 -3.7957e+03 6.5545e+03 -9.0353e+03 1.0045e+04 -9.0353e+03 6.5545e+03 -3.7957e+03 1.7213e+03 -591.4661 145.4332 -22.9073 1.7463];
%A = [1 -14.0291 95.2770 -414.6136 1.2915e+03 -3.0494e+03 5.6404e+03 -8.3317e+03 9.9295e+03 -9.5782e+03 7.4544e+03 -4.6331e+03 2.2558e+03 -832.5102 219.9203 -37.2233 3.0496];
% Pasabajos tc = tan(0.9425/2), N=13
%B = [2.9694e-09 3.8602e-08 2.3161e-07 8.4924e-07 2.1231e-06 3.8216e-06 5.0955e-06 5.0955e-06 3.8216e-06 2.1231e-06 8.4924e-07 2.3161e-07 3.8602e-08 2.9694e-09];
%A = [-1 8.7745 -36.0375 91.6024 -160.5951 204.8724 -195.4996 141.1934 - 77.1307 31.4600 -9.3083 1.8911 -0.2364 0.0137];

% Pasabajos tc = tan(0.9425/2), N=8 FUNCIONA
B = [5.6266e-06 4.5013e-05 1.5755e-04 3.1509e-04 3.9386e-04 3.1509e-04 1.5755e-04 4.5013e-05 5.6266e-06];
A = [1 -5.3910 13.0240 -18.3415 16.4225 -9.5532 3.5201 -0.7501 0.0707];

[HA, HAf] = freqz(B, A, 8192, 2*pi);

impulse_N = 81;  
impulse = [1; zeros(N-1, 1)];
hab = filter(B, A, impulse);

y_IIR = filter(B, A, x);
%audiowrite("Audios\x_IIR_filtered.wav", y_IIR, fs);

%% Graficos filtro IIR
figure(3);

subplot(4, 1, 1)
plot(0:impulse_N-1, hab);
xlabel('n (Samples)');
ylabel('Amplitude');
title('Respuesta impulsiva IIR');
grid on;

subplot(4, 1, 2)
plot(HAf, abs(HA));
title("Respuesta en frecuencia H(n): Tanto por uno")
grid on;

subplot(4, 1, 3)
plot(HAf, 20*log10(abs(HA)));
title("Respuesta en frecuencia H(n): DB")
grid on;

subplot(4, 1, 4)
plot(HAf, angle(HA));
title("Fase H(e^(j*theta))")
grid on;