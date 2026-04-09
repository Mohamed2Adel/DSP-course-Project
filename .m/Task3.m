%% ============================================================
%  Butterworth Filter Transformations
%  LPF → HPF (Coefficient Rotation)
%  LPF → BPF (Rotation by ±π/2 in Z-plane)
% ============================================================

clear; clc; close all;

%% -------------------------------------------------------------
% Load audio
%% -------------------------------------------------------------
[x, Fs] = audioread("sample_audio_file.wav");
x = x(:,1);
x = x(:);

fprintf("Audio loaded: Fs = %d Hz\n", Fs);

Nfft = 16384*2;
omega = linspace(-pi, pi, Nfft);   % digital frequency axis (rad/sample)
f_center = (-Nfft/2:Nfft/2-1) * (Fs/Nfft);

%% -------------------------------------------------------------
% LPF Specifications (Butterworth ONLY)
%% -------------------------------------------------------------
fp  = 3000;
fsb = 4000;
Ap  = 1;
As  = 50;

Wp = fp/(Fs/2);
Ws = fsb/(Fs/2);

%% -------------------------------------------------------------
% Butterworth LPF (BASE FILTER – UNCHANGED)
%% -------------------------------------------------------------
[Nb, Wnb] = buttord(Wp, Ws, Ap, As);
[b_lpf, a_lpf] = butter(Nb, Wnb, 'low');

[z_lpf, p_lpf, k_lpf] = tf2zpk(b_lpf, a_lpf);

[Hwhole_lpf, ~] = freqz(b_lpf, a_lpf, Nfft, 'whole');
Hc_lpf = fftshift(Hwhole_lpf);

%% -------------------------------------------------------------
% LPF Plots
%% -------------------------------------------------------------
figure('Name','LPF','Position',[100 50 1400 900]);

subplot(2,1,1);
zplane(b_lpf, a_lpf);
grid on;
title(sprintf("LPF – Pole-Zero (Order=%d)", Nb));

subplot(2,1,2);
plot(f_center/1e3, 20*log10(abs(Hc_lpf)+eps),'LineWidth',1.2);
grid on;
xlabel("kHz"); ylabel("dB");
title("LPF Magnitude (centered)");

%% =============================================================
% LPF to HPF Conversion (Coefficient Rotation)
%% =============================================================
n_b = 0:length(b_lpf)-1;
n_a = 0:length(a_lpf)-1;

b_hpf = b_lpf .* ((-1).^n_b);
a_hpf = a_lpf .* ((-1).^n_a);

[Hwhole_hpf, ~] = freqz(b_hpf, a_hpf, Nfft, 'whole');
Hc_hpf = fftshift(Hwhole_hpf);
gd_hpf = fftshift(grpdelay(b_hpf, a_hpf, Nfft, 'whole'));

[h_imp_hpf, n_imp] = impz(b_hpf, a_hpf, 256);

y_hpf = filter(b_hpf, a_hpf, x);
MSE_hpf = mean((x - y_hpf).^2);

%% -------------------------------------------------------------
% HPF Plots
%% -------------------------------------------------------------
figure('Name','HPF (Coefficient Rotation)','Position',[100 50 1400 900]);

subplot(3,2,1);
zplane(b_hpf, a_hpf);
grid on;
title(sprintf("HPF – Pole-Zero (Order=%d)", Nb));

subplot(3,2,2);
plot(f_center/1e3,20*log10(abs(Hc_hpf)+eps),'LineWidth',1.2);
grid on;
xlabel("kHz"); ylabel("dB");
title("HPF Magnitude (centered)");

subplot(3,2,3);
plot(f_center/1e3,unwrap(angle(Hc_hpf)),'LineWidth',1.1);
grid on;
xlabel("kHz"); ylabel("rad");
title("HPF Phase (centered)");

subplot(3,2,4);
plot(omega,gd_hpf,'m','LineWidth',1.2);
grid on;
xlabel('\omega (rad/sample)');
title("HPF Group Delay -3.14<w<3.14");

subplot(3,2,5);
Nshow = 200;
plot(n_imp(1:Nshow),h_imp_hpf(1:Nshow),'o-','LineWidth',1.2,'MarkerSize',3);
grid on;
xlabel("n"); ylabel("Amplitude");
title("HPF Impulse Response (Zoomed)");

fprintf("HPF Order = %d | MSE = %.4e\n", Nb, MSE_hpf);

%% =============================================================
% LPF to BPF Conversion (Rotation by ±π/2)
%% =============================================================
phi = pi/2;

z1 = z_lpf;     % zeros remained 
p1 = p_lpf * exp( 1j*phi); % poles rotated by pi/2

z2 = z_lpf * exp(-2j*phi); % zeros rotated by -pi 
p2 = p_lpf * exp(-1j*phi); % poles rotated by -pi/2

z_bpf = [z1 ; z2];
p_bpf = [p1 ; p2];

[b_bpf,a_bpf] = zp2tf(z_bpf,p_bpf,1);

Hnorm = polyval(b_bpf, 1j) / polyval(a_bpf, 1j);
k_bpf = 1 / abs(Hnorm);
[b_bpf,a_bpf] = zp2tf(z_bpf,p_bpf,k_bpf);

b_bpf = real(b_bpf);
a_bpf = real(a_bpf);

[Hwhole_bpf, ~] = freqz(b_bpf, a_bpf, Nfft, 'whole');
Hc_bpf = fftshift(Hwhole_bpf);
gd_bpf = fftshift(grpdelay(b_bpf, a_bpf, Nfft, 'whole'));

[h_imp_bpf, n_imp] = impz(b_bpf, a_bpf, 256);

y_bpf = filter(b_bpf, a_bpf, x);
MSE_bpf = mean((x - y_bpf).^2);

%% -------------------------------------------------------------
% BPF Plots
%% -------------------------------------------------------------
figure('Name','BPF (Rotation π/2)','Position',[100 50 1400 900]);

subplot(3,2,1);
zplane(b_bpf, a_bpf);
grid on;
title(sprintf("BPF – Pole-Zero (Order=%d)", length(a_bpf)-1));

subplot(3,2,2);
plot(f_center/1e3,20*log10(abs(Hc_bpf)+eps),'LineWidth',1.2);
grid on;
xlabel("kHz"); ylabel("dB");
title("BPF Magnitude (centered)");

subplot(3,2,3);
plot(f_center/1e3,unwrap(angle(Hc_bpf)),'LineWidth',1.1);
grid on;
xlabel("kHz"); ylabel("rad");
title("BPF Phase (centered)");

subplot(3,2,4);
plot(omega,gd_bpf,'m','LineWidth',1.2);
grid on;
xlabel('\omega (rad/sample)');
title("BPF Group Delay -3.14<w<3.14");

subplot(3,2,5);
plot(n_imp(1:Nshow),h_imp_bpf(1:Nshow),'o-','LineWidth',1.2,'MarkerSize',3);
grid on;
xlabel("n"); ylabel("Amplitude");
title("BPF Impulse Response (Zoomed)");

fprintf("BPF Order = %d | MSE = %.4e\n", length(a_bpf)-1, MSE_bpf);