clc; clear all; close all;
%% BS Filter Specs
omega_s1 = 0.240 * pi;
omega_s2 = 0.365 * pi;
omega_p2 = 0.380 * pi;
omega_p1 = 0.225* pi;
transition_bw = 0.0471;
delta = 0.15;
omega_c1 = (omega_p1 + omega_s1)/2;
omega_c2 = (omega_p2 + omega_s2)/2;
%% Kaiser window parameters
A = -20 * log10(delta);
min_width = ceil( 1 + ((A - 8)/(2.285 * transition_bw)));

alpha = -1;
if A < 21
    alpha = 0;
elseif A >=21 && A <= 50idea
    alpha = 0.5842 * (A - 21) ^ 0.4 + 0.07886 * (A - 21);
elseif A > 50
    alpha = 0.1102 * (A - 8.7);
else
    
end
beta = alpha /min_width;
%% Magnitude Response
M = min_width + 13;
w = kaiser(M,beta);
BSF_IDEAL = ideal_lpf(pi,M) - ideal_lpf(omega_c2,M) + ideal_lpf(omega_c1,M);
BSF_FIR = BSF_IDEAL .* w';
[H,f] = freqz(BSF_FIR,1,1024, 400e3);
figure();
plot(f, abs(H));
hold on;
xline(48e3, 'm--', 'LineWidth', 1.5);
hold on;
xline(73e3, 'm--', 'LineWidth', 1.5);
hold on;
xline(45e3, 'm--', 'LineWidth', 1.5);
hold on;
xline(76e3, 'm--', 'LineWidth', 1.5);
hold on;
yline(1.15, 'r--', 'LineWidth', 1.5);
hold on;
yline(0.85, 'r--', 'LineWidth', 1.5);
hold on;
yline(0.15, 'r--', 'LineWidth', 1.5); 
xlabel('f in 10^4 Hz');
ylabel('|H(e^{j2 \pi f})');
title('Magnitude Response of Discrete Time FIR Bandstop Filter');
set(gca, 'XTick', [45e3, 48e3, 73e3, 76e3], 'xticklabel', {'f_{p1}', 'f_{s1}', 'f_{s2}', 'f_{p2}'});
set(gca, 'YTick', [0.15, 0.85, 1, 1.15], 'yticklabel', {'\delta_2 = 0.15', '1 - \delta_1 = 0.85', '1', '1 + \delta_1 = 1.15'});
%% phase and impulse response 
disp(BSF_FIR);
fvtool(BSF_FIR, 'Analysis', 'Phase');
fvtool(BSF_FIR, 'Analysis', 'Impulse');
%% ideal lpf code
function lpass = ideal_lpf(wc,M);

alpha = (M-1)/2;
n = [0:1:(M-1)];
m = n - alpha + eps;
lpass = sin(wc*m) ./ (pi*m);
end