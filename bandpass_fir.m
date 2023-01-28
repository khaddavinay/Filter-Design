clc; clear all; close all;
%% BP Filter Specifications
omega_s1 = 0.2518 * pi;
omega_s2 = 0.4407 * pi;
omega_p2 = 0.4296 * pi;
omega_p1 = 0.2629 * pi;
transition_bw = 0.0111 * pi;
delta = 0.15;
omega_c1 = (omega_p1 + omega_s1)/2;
omega_c2 = (omega_p2 + omega_s2)/2;
%% Kaiser window parameters
A = -20 * log10(delta);
min_width = ceil( 1 + ((A - 8)/(2.285 * transition_bw)));

alpha = -1;
if A < 21
    alpha = 0;
elseif A >=21 && A <= 50
    alpha = 0.5842 * (A - 21) ^ 0.4 + 0.07886 * (A - 21);
elseif A > 50
    alpha = 0.1102 * (A - 8.7);
else
    
end
beta = alpha /min_width;
%% kaiser window multiplication
M = min_width + 3;
K = kaiser(M,beta);
BPF_IDEAL = ideal_lpf(omega_c2, M) - ideal_lpf(omega_c1, M);
BPF_FIR= BPF_IDEAL .* K';
[H,f] = freqz(BPF_FIR,1,1024, 540e3);
figure();
plot(f, abs(H));
hold on;
xline(71e3, 'm--', 'LineWidth', 1.5);
hold on;
xline(116e3, 'm--', 'LineWidth', 1.5);
hold on;
xline(68e3, 'm--', 'LineWidth', 1.5);
hold on;
xline(119e3, 'm--', 'LineWidth', 1.5);
hold on;
yline(1.15, 'r--', 'LineWidth', 1.5);
hold on;
yline(0.85, 'r--', 'LineWidth', 1.5);
hold on;
yline(0.15, 'r--', 'LineWidth', 1.5); 
xlabel('f in 10^4 Hz');
ylabel('|H(e^{j2 \pi f})');
title('Magnitude Response of Discrete Time FIR Bandpass Filter');
set(gca, 'XTick', [68e3, 71e3, 116e3, 119e3], 'xticklabel', {'f_{s1}', 'f_{p1}', 'f_{p2}', 'f_{s2}'});
set(gca, 'YTick', [0.15, 0.85, 1, 1.15], 'yticklabel', {'\delta_2 = 0.15', '1 - \delta_1 = 0.85', '1', '1 + \delta_1 = 1.15'});
%% displays the impulse and the phase response 
fvtool(BPF_FIR, 'Analysis', 'Phase');
fvtool(BPF_FIR, 'Analysis', 'Impulse');

%% ideal lpf code
function lpass = ideal_lpf(wc,M);

alpha = (M-1)/2;
n = [0:1:(M-1)];
m = n - alpha + eps;
lpass = sin(wc*m) ./ (pi*m);
end