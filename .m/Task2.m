%% ============================================================
%  Task 2 – Digital IIR LPF Design & Analysis (Centered Plots)
%  Filter types: Butterworth, Chebyshev I, Chebyshev II, Elliptic
%  All frequency plots from -Fs/2 to +Fs/2
% ============================================================

clear; clc; close all;

%% -------------------------------------------------------------
% Load audio
%% -------------------------------------------------------------
[x, Fs] = audioread("sample_audio_file.wav");
x = mean(x,2);  % averaging 2 channels like in task 1
%x = x(:);

fprintf("Audio loaded: Fs = %d Hz\n", Fs);

%% -------------------------------------------------------------
% Filter specifications
%% -------------------------------------------------------------
fp  = 3000;    % Passband edge
fsb = 4000;    % Stopband edge
Ap  = 1;       % Passband ripple (dB)
As  = 50;      % Stopband attenuation (dB)

Wp = fp/(Fs/2);
Ws = fsb/(Fs/2);

%% -------------------------------------------------------------
% Minimum-order filter designs
%% -------------------------------------------------------------
[Nb, Wnb]       = buttord(Wp, Ws, Ap, As);
[b_butt, a_butt] = butter(Nb, Wnb, "low");

[Nc1, Wnc1]     = cheb1ord(Wp, Ws, Ap, As);
[b_ch1, a_ch1]   = cheby1(Nc1, Ap, Wnc1, "low");

[Nc2, Wnc2]     = cheb2ord(Wp, Ws, Ap, As);
[b_ch2, a_ch2]   = cheby2(Nc2, As, Wnc2, "low");

[Ne, Wne]       = ellipord(Wp, Ws, Ap, As);
[b_ell, a_ell]   = ellip(Ne, Ap, As, Wne, "low");

%% -------------------------------------------------------------
% Filter the audio
%% -------------------------------------------------------------
y_butt = filter(b_butt, a_butt, x);
y_ch1  = filter(b_ch1 , a_ch1 , x);
y_ch2  = filter(b_ch2 , a_ch2 , x);
y_ell  = filter(b_ell , a_ell , x);

%% -------------------------------------------------------------
% Prepare filter info for looping
%% -------------------------------------------------------------
filters = {
    "Butterworth", b_butt, a_butt, y_butt, Nb;
    "Chebyshev I", b_ch1,  a_ch1,  y_ch1,  Nc1;
    "Chebyshev II",b_ch2,  a_ch2,  y_ch2,  Nc2;
    "Elliptic",    b_ell,  a_ell,  y_ell,  Ne
};

Nfft = 32768;         % FFT points for dense frequency response
results = struct();

%% =============================================================
% Loop over filters
%% =============================================================
for k = 1:size(filters,1)

    name  = filters{k,1};
    b     = filters{k,2};
    a     = filters{k,3};
    y     = filters{k,4};
    order = filters{k,5};

    fprintf("\n=== %s Filter ===\n", name);
    fprintf("Order = %d\n", order);

    %% ------------------- Frequency Response (centered) -------------------
    [H_whole, ~] = freqz(b,a,Nfft,'whole');  % 0..2π
    Hs = fftshift(H_whole);                  % center → -π..π

    f_center = (-Nfft/2:Nfft/2-1)*(Fs/Nfft); % -Fs/2 → Fs/2 (Hz)

    %% ------------------- Phase -------------------
    phase = angle(Hs);
    phase_un = unwrap(phase); % what unwrap do it plot every pi take a magnitude and then connect the point together

    %% ------------------- Group Delay (centered -pi to pi) -------------------
    [gd, ~] = grpdelay(b,a,Nfft,'whole');  % 0..2π
    gd = fftshift(gd);                     % shift to -π..π

    %% ------------------- Impulse Response -------------------
    imp_len = 256;
    [h_imp, n_imp] = impz(b,a,imp_len);

    %% ------------------- Metrics -------------------
    mse = mean((x - y).^2);
    Ein = sum(x.^2);
    Eout = sum(y.^2);
    loss = (1 - Eout/Ein)*100;
    if loss < 0, loss = 0; end

    results(k).name = name;
    results(k).order = order;
    results(k).mse = mse;
    results(k).loss = loss;

    fprintf("MSE = %.4e | Energy Lost = %.4f %%\n", mse, loss);

    %% =========================================================
    %     ONE FIGURE PER FILTER (all subplots inside)
    %% =========================================================
    figure('Name',name,'NumberTitle','off','Position',[100 50 1400 900]);

    %-----------------------------------------------------------
    % 1) Pole-Zero
    %-----------------------------------------------------------
    subplot(3,2,1);
    zplane(b,a);
    title([name " – Pole-Zero"]);
    grid on;

    %-----------------------------------------------------------
    % 2) Magnitude Response (centered)
    %-----------------------------------------------------------
    subplot(3,2,2);
    plot(f_center/1e3, 20*log10(abs(Hs)+eps),'LineWidth',1.2); 
    grid on;
    xlabel("Frequency (kHz)");
    ylabel("Magnitude (dB)");
    title("Magnitude Response (centered)");

    %-----------------------------------------------------------
    % 3) Phase Response (centered, unwrapped)
    %-----------------------------------------------------------
    subplot(3,2,3);
    plot(f_center/1e3, phase_un,'LineWidth',1.1);
    grid on;
    xlabel("Frequency (kHz)");
    ylabel("Phase (rad)");
    title("Phase Response (centered)");

    %-----------------------------------------------------------
    % 4) Group Delay (centered, -π ≤ w ≤ π)
    %-----------------------------------------------------------
    subplot(3,2,4);
    omega = linspace(-pi, pi, Nfft);   % rad/sample
    plot(omega, gd, 'm','LineWidth',1.2);
    grid on;
    xlabel("\omega (rad/sample)");
    ylabel("Group Delay (samples)");
    title("Group Delay ( -3.14<w<3.14 )");
    

    %-----------------------------------------------------------
    % 5) Impulse Response
    %-----------------------------------------------------------
    subplot(3,2,5);
    stem(n_imp, h_imp,'filled');
    grid on;
    xlabel("n");
    ylabel("h[n]");
    title("Impulse Response");

end

%% =============================================================
% PRINT SUMMARY
%% =============================================================
fprintf("\n==================== SUMMARY ====================\n");
for k = 1:4
    fprintf("%s: Order=%d | MSE=%.3e | Lost=%.2f%%\n", ...
        results(k).name, results(k).order, results(k).mse, results(k).loss);
end

[~,idx1] = min([results.mse]);
[~,idx2] = min([results.loss]);

fprintf("\nBest (MSE): %s\n", results(idx1).name);
fprintf("Best (Energy): %s\n", results(idx2).name);
fprintf("=================================================\n\n");
