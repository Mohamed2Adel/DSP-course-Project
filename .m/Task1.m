%                                                 %%% TASK 1 %%%
%                                                 % H(z)=1+0.9z^(-1000)+0.8z^(-2000)+0.7z^(-3000)

clc; clear; close all;

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
[gd, w_gd] = grpdelay(b,a,2048,'whole');       % w_gd: 0 → 2π
w_gd_shift = w_gd - pi;                       % shift to -π → π

figure;
plot(w_gd_shift, gd, 'LineWidth', 1.1);
xlabel('Frequency (rad/sample)'); ylabel('Group Delay (samples)');
title('Group Delay');
grid on;

%% ==========================================================
% 5. IMPULSE RESPONSE
%% ==========================================================
N = 4000;
h = impz(b,a,N);

figure;
stem(h, 'filled');
xlabel('n'); ylabel('h[n]');
title('Impulse Response');
grid on;

