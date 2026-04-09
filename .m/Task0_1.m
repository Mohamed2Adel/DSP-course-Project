clc; clear; close all;
                                             %% 
                                                %%%%%%%% Task 0 %%%%%%%%

%%% plot time and frequency domain of the Audio file and finding the energy from both domains %%%%
%% 1. Read Audio File
[x_n, Fs] = audioread('sample_audio_file.wav');   % x: audio signal, Fs: sampling frequency
x_n = mean(x_n, 2);  % convert a stereo into mono signal so we can do conv of vector with vector

fprintf('Sampling frequency Fs = %.2f Hz\n', Fs);
N = length(x_n);
t = (0:N-1)/Fs;

%% 2. Plot time-domain signal x[n]
figure;
plot(t, x_n);
xlabel('Time (sec)'); ylabel('Amplitude');
title('Time-Domain Signal x[n]');
grid on;

%% 3. Plot magnitude spectrum |X(f)|
X = fftshift(fft(x_n));
f = (-N/2:N/2-1) * (Fs/N);

figure;
plot(f, abs(X));
xlabel('Frequency (Hz)'); ylabel('|X(f)|');
title('Magnitude Spectrum of x[n]');
grid on;

%% 4. Compute total energy (time domain)
E_time = sum(abs(x_n).^2);

%% 5. Compute total energy (frequency domain)
% Parseval: E = (1/N) * sum(|X[k]|^2)
E_freq = (1/N) * sum(abs(X).^2);

fprintf('\nTotal energy in time domain  = %.4f\n', E_time);
fprintf('Total energy in freq domain  = %.4f\n', E_freq);
fprintf('Difference (should be very small) = %.4e\n', abs(E_time - E_freq));








%                                                 %%% TASK 1 %%%
%                                                 % H(z)=1+0.9z^(-1000)+0.8z^(-2000)+0.7z^(-3000)


%% Define filter H(z)
Fs = 48000;
b = [1 zeros(1,999) 0.9 zeros(1,999) 0.8 zeros(1,999) 0.7];
a = 1;

%% Frequency response (whole range, then shift to -Fs/2 to Fs/2)
[H,w] = freqz(b,a,4096,'whole');        % w: 0 → 2π
f = (w - pi) * Fs / (2*pi);             % convert to -Fs/2 → Fs/2
H_shift = fftshift(H);                  % shift response

%% Magnitude in dB
mag_dB = 20*log10(abs(H_shift));

%% Phase in dB (requested)
phase = angle(H_shift);
phase_dB = 20*log10(abs(phase) + 1e-12);

%% ==========================================================
% 1. POLE-ZERO PLOT
%% ==========================================================
figure;
zplane(b,a);
title('Pole-Zero Plot');

%% ==========================================================
% 2. MAGNITUDE (dB) + PHASE (dB)  over (-Fs/2 , Fs/2)
%% ==========================================================
figure;

subplot(2,1,1);
plot(f, mag_dB, 'LineWidth', 1.1);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Magnitude Response (dB)');
grid on;

subplot(2,1,2);
plot(f, phase_dB, 'LineWidth', 1.1);
xlabel('Frequency (Hz)'); ylabel('Phase (dB)');
title('Phase Response (dB)');
grid on;

%% ==========================================================
% 3. PHASE RESPONSE ONLY (RADIANS) over (-Fs/2 , Fs/2)
%% ==========================================================
figure;
plot(f, unwrap(angle(H_shift)), 'LineWidth', 1.1);
xlabel('Frequency (Hz)'); ylabel('Phase (rad)');
title('Phase Response (Radians)');
grid on;

%% ==========================================================
% 4. GROUP DELAY over (-π , π)
%% ==========================================================
[gd, w_gd] = grpdelay(b,a,4096,'whole');       % w_gd: 0 → 2π
w_gd_shift = w_gd - pi;                       % shift to -π → π

figure;
plot(w_gd_shift, gd, 'LineWidth', 1.1);
xlabel('Frequency (rad/sample)'); ylabel('Group Delay (samples)');
title('Group Delay');
grid on;

%% ==========================================================
% 5. IMPULSE RESPONSE
%% ==========================================================
N = 14000;
h = impz(b,a,N);

figure;
stem(h, 'filled');
xlabel('n'); ylabel('h[n]');
title('Impulse Response');
grid on;



%% ==========================================================
% Apply digital echo system using convolution
%% ==========================================================
% h(n) is already defined from impz(b,a,N) or you can use:
% h = b;  % FIR filter coefficients themselves are the impulse response

y_n = conv(x_n, b);   % convolution of input with impulse response h=b

%% Compute MSE
N_x = length(x_n);                  % length of original signal
% Compare only first N_x samples since convolution output is longer
MSEy1 = (1/(N_x+1)) * sum( (y_n(1:N_x) - x_n).^2 );

%% Display result
fprintf('Mean-Square Error (MSE) of y_n: %.6f\n', MSEy1);

%% ==========================================================
%                %%% TASK 1 — Equalizer G(z) = 1/H(z) %%%
% ===========================================================

% H(z) coefficients from previous section
bH = b;   % numerator of H(z)
aH = a;   % denominator (1)

%% Define G(z) = 1 / H(z)
bG = aH;      % numerator = 1
aG = bH;      % denominator = H(z) coefficients

%% Frequency response of G(z)
[HG, wG] = freqz(bG, aG, 4096, 'whole');    % wG from 0 to 2pi
fG = (wG - pi) * Fs / (2*pi);               % shift it from -Fs/2 to Fs/2
HG_shift = fftshift(HG);                    % shift response 

%% Magnitude (dB)
magG_dB = 20*log10(abs(HG_shift));

%% Phase (dB)
phaseG = angle(HG_shift);
phaseG_dB = 20*log10(abs(phaseG) + 1e-12);

%% ==========================================================
% 1. POLE-ZERO PLOT of G(z)
%% ==========================================================
figure;
zplane(bG, aG);
title('G(z) — Pole-Zero Plot (Inverse Filter)');

%% ==========================================================
% 2. MAGNITUDE + PHASE in dB
%% ==========================================================
figure;

subplot(2,1,1);
plot(fG, magG_dB, 'LineWidth', 1.1);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('G(z) Magnitude Response (dB)');
grid on;

subplot(2,1,2);
plot(fG, phaseG_dB, 'LineWidth', 1.1);
xlabel('Frequency (Hz)'); ylabel('Phase (dB)');
title('G(z) Phase Response (dB)');
grid on;

%% ==========================================================
% 3. PHASE RESPONSE (Radians)
%% ==========================================================
figure;
plot(fG, unwrap(angle(HG_shift)), 'LineWidth', 1.1);
xlabel('Frequency (Hz)'); ylabel('Phase (rad)');
title('G(z) Phase Response (Radians)');
grid on;

%% ==========================================================
% 4. GROUP DELAY
%% ==========================================================
[gdG, wGD] = grpdelay(bG, aG, 2048, 'whole');
wGD_shift = wGD - pi;

figure;
plot(wGD_shift, gdG, 'LineWidth', 1.1);
xlabel('Frequency (rad/sample)'); ylabel('Group Delay (samples)');
title('G(z) Group Delay');
grid on;

%% ==========================================================
% 5. IMPULSE RESPONSE of G(z)
%% ==========================================================
NG = 40000;
g_n = impz(bG, aG, NG);

figure;
stem(g_n, 'filled');
xlabel('n'); ylabel('g[n]');
title('G(z) Impulse Response');
grid on;


%% Apply equalizer system G to the echoed signal y_n
y2_n = conv(y_n, g_n);      % Output of the equalizer

%% Compute MSE between y2_n and original x_n
% Make the lengths match (use the first N_x samples)
MSEy2 = (1/(N_x+1)) * sum( (y2_n(1:N_x) - x_n).^2 );

%% Display result
fprintf('Mean-Square Error (MSE) of y2_n: %.6f\n', MSEy2);






