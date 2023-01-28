clc; close all;
%% FILTER SPECIFICATIONS
wp1 = 0.706; 
ws1 = 0.754; 
ws2 = 1.146; 
wp2 = 1.194;
op1= 0.3689194771098271; 
os1 = 0.3959280087977212;
os2 = 0.645691821042554;
op2 = 0.6795992982245265;
%% Pole plot of LPF
N = 5;
Omega_0 = 0.5007;
B = 0.31067;
%% EQUIVALENT LOWPASS FILTER SPECIFICATIONS
D2 = 43.44;
D1 = 0.384;
Omega_ls = 1.2069;
Omega_lp = 1;
epsilon = sqrt(D1);
N = ceil((acosh(sqrt(D2/D1))) / (acosh(Omega_ls / Omega_lp)));
fprintf('N = %d\n', N);
%% LHP POLES 
a = 0:2 * N - 1;
A = (2 * a + 1) .* (pi /(2 * N));
B = (1 / N) * asinh(1 / epsilon);
poles = (-1 * sin(A) * sinh(B)) + 1i * (cos(A) * cosh(B));
figure();
plot(poles, '.', 'MarkerSize', 20);
daspect([1 1 1]);
hold on;
x = linspace(-pi, pi, 10000);
a = sinh(B) .* cos(x);
b = cosh(B) .* sin(x);
plot(a,b);
hold on;
plot(0,0,'g*');
LHP = [];
for i=1:size(poles,2)
    hold on;
    plot([0, real(poles(1,i))], [0, imag(poles(1,i))], 'r-');
    if real(poles(1,i)) < 0
        disp(poles(1,i));
        LHP = [LHP, poles(1,i)];
    end
end
plot([0, 0],[-1.5, 1.5], 'k-');
plot([-1.5, 1.5],[0, 0], 'k-');
title('Poles of $H_{LPF}(s_L) \, H_{LPF}(-s_L)$ in the $s_L$ plane', 'Interpreter', 'latex');
xlabel('\Sigma_k');
ylabel('\Omega_k');
%% ANALOG LOWPASS TRansfer function
syms s;
den=1;
num =1;
for i=1:size(LHP, 2)
    num = num * LHP(1,i);
    den = den * (s - LHP(1,i));
end
H_LPF = num/den;
num_f = num;
den_f = den;
den = sym2poly(den);
w = linspace(0,2.5,1000);
h = freqs(num,den,w);
figure();
plot(w,abs(h),'LineWidth',1);
hold on;
yline(0.85,['--','r']);
xline(1,['--','g']);
xline(1.2069,['--','g']);
yline(1.,['--','r']);
yline(0.15,['--','r']);
ylim([0,1.2]);
title("Equivalent Chebyshew Low Pass Magnitude Response");
xlabel("\Omega_L");
ylabel("|H(j\Omega_L|");
set(gca,'Fontsize',5,'XTick', [0, 1, 1.2069], 'xticklabel', {'0', '\Omega_{Lp} = 1','\Omega_{Ls} = 1.2069'});
set(gca, 'YTick', [0, 0.15, 0.85,1], 'yticklabel', {'0', '\delta = 1','1-\delta = 0.85','1'});
figure();
plot(w,angle(h));
title("Equivalent Chebyshew Low Pass Phase Response");
xlabel("\Omega_L");
ylabel("\angle H(j\Omega_L");

%% ANALOG BANDPASS TRANSFORMATION

syms s1;
Omega_0 = 0.5007;
B = 0.31067;
s1 = (B * s)./(s^2 + Omega_0^2);
h = subs(den_f, s, s1);
H_BSF = subs(H_LPF,s,s1);
H = num_f/h;

%% MAGNITUDE AND PHASE PLOT
[num,den] = numden(H);
k = subs(num(1),s,1);
num_coeff = sym2poly(num/k);
den_coeff = sym2poly(den/k);
% disp('************************************************');
% disp(num_coeff);
% disp(den_coeff);
disp('************************************************');
w = linspace(0,2,1000);
[g,w] = freqs(num_coeff,den_coeff,w);
%disp(g);
figure();
plot(w,abs(g),'LineWidth',1);
hold on;
yline(0.85,['--','r']);
xline(os1,['--','g']);
xline(os2,['--','g']);
xline(op1,['--','g']);
xline(op2,['--','g']);
yline(1.,['--','r']);
yline(0.15,['--','r']);
ylim([0,1.2]);
xlim([0,1])
title("Chebyshew Bandstop Filter Magnitude Response");
xlabel("\Omega");
ylabel("|H(j\Omega|");
set(gca, 'XTick', [op1, os1, os2, op2], 'xticklabel', {'\Omega_{p1}', '\Omega_{s1}', '\Omega_{s2}', '\Omega_{p2}'});
set(gca, 'YTick', [0.15, 0.85, 1], 'yticklabel', {'\delta_2 = 0.15', '1 - \delta_1 = 0.85', '1'});
%% phase response
figure();
plot(w,angle(g),'LineWidth',1);
hold on;
title("Chebyshew Bandstop Filter Phase Response");
xlabel("\Omega");
ylabel('\angle H(j \Omega)');
%% Discrete domain transformation
syms z;
c = subs(H_BSF, s, (z-1) / (z+1));
disp(c);

[num,den] = numden(c);
num_coeff = sym2poly(num);
den_coeff = sym2poly(den);
% disp('************************************************');
disp(num_coeff);
disp(den_coeff);
w = linspace(0,2,1024);
[g,w] = freqz(num_coeff,den_coeff,1024*1024,400e3);
figure();
plot(w,abs(g),'Linewidth',1);
hold on;
yline(0.85,['--','r']);
xline(45e3,['--','g']);
xline(48e3,['--','g']);
xline(73e3,['--','g']);
xline(76e3,['--','g']);
yline(1.,['--','r']);
yline(0.15,['--','r']);
ylim([0,1.2]);
xlim([0,140e3])
title("Discrete Time Bandstop Filter Magnitude Response");
xlabel("\Omega");
ylabel("|H(j\Omega|");
set(gca, 'XTick', [45e3, 48e3, 73e3, 76e3], 'xticklabel', {'f_{p1}', 'f_{s1}', 'f_{s2}', 'f_{p2}'});
set(gca, 'YTick', [0.15, 0.85, 1, 1.15], 'yticklabel', {'\delta_2 = 0.15', '1 - \delta_1 = 0.85', '1', '1 + \delta_1 = 1.15'});
figure();
plot(w,angle(g),'LineWidth',1);
hold on;
title("Discrete Time Bandstop Filter Phase Response");
xlabel("\Omega");
ylabel('\angle H(j \Omega)');
















