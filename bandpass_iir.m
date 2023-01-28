clc; close all;
tic;
%% Poles of H_{analog}(s_L) H_{analog}(-s_L)
N = 21;
Omega_c = 1.02462;
a = 1:2*N;
poles = Omega_c .* exp(1i .* (pi/2) .* (1 + (2.*a + 1)./ N));
figure();
plot(poles, '*', 'MarkerSize', 10);
hold on;
x = linspace(-pi, pi, 10000);
a = Omega_c .* cos(x);
b = Omega_c .* sin(x);
plot(a,b);
hold on;
plot(0,0,'r*');
LHP = [];
for i=1:size(poles,2)
    hold on;
    plot([0, real(poles(1,i))], [0, imag(poles(1,i))], 'r-');
    if real(poles(1,i)) < 0
%         disp(poles(1,i));
        LHP = [LHP, poles(1,i)];
    end
end
plot([0, 0],[-1.5, 1.5], 'k-');
plot([-1.5, 1.5],[0, 0], 'k-');
title('Poles of $H_{LPF}(s_L) \, H_{LPF}(-s_L)$ in the $s_L$ plane', 'Interpreter', 'latex');
xlabel('\Sigma_k');
ylabel('\Omega_k');
%% Finding poles of H_{analog}(s_L)
syms s x;
temp = 1;
for i=1:size(LHP,2)
    temp = temp * (s - LHP(1,i));
end
%% Magnitude and Phase response of H_{analog}(s_L) 
den_coeff = sym2poly(temp);
num = 1.665;
w = linspace(0,2,1000);
h = freqs(num,den_coeff,w);
figure();
plot(w,abs(h));
hold on;
fplot(s - s - 0.15 + 1, 'k-', 'Markersize', 10);
hold on;
fplot(x-x+1, 'k-', 'Markersize', 10);
hold on;
fplot(x - x + 0.15, 'r-', 'Markersize', 10);
xline(1, 'm-');
hold on;
xline(1.1226,'m-');
axis([0 2 0 1.2]);
set(gca,'Fontsize',5,'XTick',[0, 1, 1.1226], 'xticklabel', {'0', '\Omega_{Lp} = 1','\Omega_{Ls} = 1.1226'});
set(gca, 'YTick', [0.15, 0.85, 1], 'yticklabel', {'\delta_2 = 0.15', '1 - \delta_1 = 0.85', '1'});
daspect([1 1 1]);
title('Equivalent Butterworth Lowpass filter mangitude response');
xlabel('\Omega_L');
ylabel('|H(j \Omega_L)|');
figure();
plot(w,angle(h));
xlabel('\Omega_L');
ylabel('\angle H(j \Omega_L)');
title('Equivalent Butterworth Lowpass filter phase response');




%% Bandpass analog frequency transformation
syms s1;
Omega_0 = 0.592;
B = 0.362;
s1 = (s^2 + Omega_0^2) ./ (B * s);
h = subs(temp, s, s1);
H = 1.6665/h;
%% Normalization of the coefficients 
[num,den] = numden(H);
k = subs(num(1),s,1);   %  numerator coefficient
num_coeff = sym2poly(num/k);
den_coeff = sym2poly(den/k);
%% magnitude and phase plot of BPF
w = linspace(0,2,1000);
[g,w] = freqs(num_coeff,den_coeff,w);
figure();
plot(w,abs(g));
set(gca,'Fontsize',4,'XTick', [0.4175, 0.4381, 0.80012, 0.8291], 'xticklabel', {'\Omega_{s1}', '\Omega_{p1}', '\Omega_{p2}', '\Omega_{s2}'});
set(gca, 'YTick', [0.15, 0.85, 1], 'yticklabel', {'\delta_2 = 0.15', '1 - \delta_1 = 0.85', '1'});
hold on;
xline(0.4175, 'm-');hold on;
xline(0.4381, 'm-');hold on;
xline(0.80012, 'm-');hold on;
xline( 0.8291, 'm-');
hold on;
fplot(s - s - 0.15 + 1, 'k-', 'Markersize', 10);
hold on;
fplot(x-x+1, 'k-', 'Markersize', 10);
hold on;
fplot(x - x + 0.15, 'r-', 'Markersize', 10);
daspect([1 1 1]);
title('Butterworth Bandpass filter mangitude response');
xlabel('\Omega');
ylabel('|H(j \Omega)|');
axis([0 1.2 0 1.2]);

figure();
hold on;
plot(w,angle(g));
title('Butterworth Bandpass filter phase response');
xlabel('\Omega');
ylabel('\angle H(j \Omega)');

%% Discrete domain transformation
syms w;
c = subs(H, s, (1-1/w) / (1+1/w));
disp(c);
[num,den] = numden(c);
k = -subs(num(1),w,0);
num_coeff = sym2poly(num/k);
den_coeff = sym2poly(den/k);
% disp(num_coeff);
% disp(den_coeff);
w = linspace(0,2,1024);
[g,w] = freqz(num_coeff,den_coeff,1024,540e3);
figure();
plot(w,abs(g));
figure();
plot(w,angle(g));