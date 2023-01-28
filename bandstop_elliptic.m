% Please execute section wise
clc; clear ; close all;
%% elliptical filter parameter specifications calculations
ws =1.2069;
wp=1;
Gp = 0.85;
D1 = 0.384;
D2 = 43.44;
Ap = 1.4116;
As = 16.478;
ep = sqrt(10^(Ap/10) - 1); % ripple factors
es = sqrt(10^(As/10) - 1);
k = wp/ws;
k1 = ep/es;
kf = sqrt(1-k^2);  
t = (k1)^2;
kf1 = sqrt(1-t);
%% JACOBIAN INTEGRAL CALCULATION
syms x;
f = 1/sqrt(1-k^2*sin(x)^2);
g = 1/sqrt(1-kf^2*sin(x)^2);
h = 1/sqrt(1-k1^2*sin(x)^2);
i = 1/sqrt(1-kf1^2*sin(x)^2);
K = double(int(f,0,pi/2));
Kf = double(int(g,0,pi/2)); 
K1 = double(int(h,0,pi/2));
Kf1 = double(int(i,0,pi/2));
disp(K);   %Jacobian integral output  
disp(Kf);
disp(K1);
disp(Kf1);
N = ceil((Kf1*K)/(K1*Kf));  %filter order 
disp(N);
%% updation of k value
u1 = 1/N;
u2= Kf1/4;
m = (Kf1*3)/4;
%kj = (((kf1)^N)*(sne(u1,kf1)^4)*(sne(u2,kf1)^4));
kj = ellipdeg(N,k1);
%% pole zero calculation for the above specifications
L = floor(N/2); r = mod(N,2); % L is the number of second-order sections
i = (1:L)'; u = (2*i-1)/N; zeta_i = cde(u,kj); 
z = wp*1j./(kj*zeta_i); % zeros of elliptic rational function
ZEROES = z;
v0 = -1j*asne(1j/ep, k1)/N; 
p = wp*1j*cde((u-1j*v0), kj); % filter poles
poles = p;
p0 = wp*1j*sne(1j*v0, kj); % first-order pole, needed when N is odd
B = [ones(L,1), -2*real(1./z), abs(1./z).^2]; % second-order numerator coefficients
A = [ones(L,1), -2*real(1./p), abs(1./p).^2]; % second-order denominator coefficients
if r==0 % prepend first-order sections
    B = [Gp, 0, 0; B]; A = [1, 0, 0; A];
else
    B = [1, 0, 0; B]; A = [1, -real(1/p0), 0; A];
end
ZEROES = cplxpair([ZEROES; conj(ZEROES)]); % append conjugate zeros
p = cplxpair([p; conj(p)]); % append conjugate poles
if r==1, p = [p; p0]; end % append first-order pole when N is odd
H0 = Gp^(1-r); % dc gain
%% Elliptic Lowpass filter transfer function and plots
syms sl Omega_L;
H_LPF = (1 + 0.6295 * sl^2)/ ((1 + 1.6047 * sl) * (1 + 0.2305 * sl + 0.9994 * sl^2));
H_LPF_freq = subs(H_LPF, sl, 1i * Omega_L);
[ns, ds] = numden(H_LPF);
nsl = sym2poly(ns);
dsl = sym2poly(ds);
k = ds(1);
nsl = nsl / k;
ds = ds / k;
disp(nsl);
disp(dsl);
figure();
fplot(abs(H_LPF_freq),'LineWidth',1);
hold on;
fplot(sl - sl - 0.15 + 1, 'g--', 'Markersize', 10);
hold on;
fplot(sl-sl+1, 'g--', 'Markersize', 10);
hold on;
fplot(sl - sl + 0.15, 'g--', 'Markersize', 10);
xline(1, 'magenta--','LineWidth',1);
hold on;
xline(1.2069,'magenta--','LineWidth',1);
axis([0 2 0 1.2]);
set(gca,'FontSize',5, 'XTick', [0, 1, 1.2069], 'xticklabel', {'0', '\Omega_{Lp} = 1','\Omega_{Ls} = 1.2069'});
set(gca, 'YTick', [0.15, 0.85, 1], 'yticklabel', {'\delta_2 = 0.15', '1 - \delta_1 = 0.85', '1'});
daspect([1 1 1]);
title('Equivalent Elliptic Lowpass filter mangitude response');
xlabel('\Omega_L');
ylabel('|H(j \Omega_L)|');
figure();
fplot(angle(H_LPF_freq));
xlabel('\Omega_L');
ylabel('\angle H(j \Omega_L)');
title('Equivalent elliptic Lowpass filter phase response');

 %% Analog lowpass to bandpass frequency transformation
 syms s Omega;
 Omega_0 = 0.5007;
 B1 = 0.31067;
 H_BSF = subs(H_LPF, sl,(B1 * s)/(s^2 + Omega_0^2));
[ns, ds] = numden(H_BSF);
ns = sym2poly(ns);
ds = sym2poly(ds);
k = ds(1);
ns = ns ./ k;
ds = ds ./ k;
disp(ns);
disp(ds);
sys = tf(ns,ds);   % Analog BSF transfer function 
H_BSF_freq = subs(H_BSF, s, 1i * Omega);
figure();
fplot(abs(H_BSF_freq),'LineWidth',1);
set(gca, 'XTick', [0.3689, 0.3959, 0.6457, 0.6795], 'xticklabel', {'\Omega_{p1}', '\Omega_{s1}', '\Omega_{s2}', '\Omega_{p2}'});
set(gca, 'YTick', [0.15, 0.85, 1], 'yticklabel', {'\delta_2 = 0.15', '1 - \delta_1 = 0.85', '1'});
hold on;
xline(0.3689, 'r--');hold on;
xline(0.3959, 'r--');hold on;
xline(0.6457, 'r--');hold on;
xline( 0.6795, 'r--');
hold on;
fplot(s - s - 0.15 + 1, 'g--', 'Markersize', 10);
hold on;
fplot(s-s+1, 'g--', 'Markersize', 10);
hold on;
fplot(s - s + 0.15, 'g--', 'Markersize', 10);
daspect([1 1 1]);
title('Elliptic Bandstop filter mangitude response');
xlabel('\Omega');
ylabel('|H(j \Omega)|');
axis([0 1.5 0 1.2]);
%-------------------------------------------------------------
%% Analog to z bilinear transformation
syms z;
Hz = subs(H_BSF, s, (z-1)/(z+1));
[Nz, Dz] = numden(Hz);
Nz = sym2poly(Nz);
Dz = sym2poly(Dz);
k = Dz(1);
Nz = Nz / k;
Dz = Dz / k;
disp(Nz);
disp(Dz);
[H,f] = freqz(Nz,Dz,1024*1024, 400e3);
SYS = tf(Nz,Dz);
figure();
plot(f,abs(H),'LineWidth',1);
axis([20e3 100e3 0 1.3]);
hold on;
yline(1.00, 'g--', 'LineWidth', 1);
hold on;
yline(0.85, 'g--', 'LineWidth', 1);
hold on;
yline(0.15, 'g--', 'LineWidth', 1);
hold on;
xline(45e3, 'r--', 'LineWidth', 1);
hold on;
xline(48e3, 'r--', 'LineWidth', 1);
hold on;
xline(73e3, 'r--', 'LineWidth', 1);
hold on;
xline(76e3, 'r--', 'LineWidth', 1);
xlabel('f in 10^4 Hz');
ylabel('|H(e^{j 2\pi f})|');
title('Magnitude Response of the Discrete Time Bandstop Filter');
set(gca, 'XTick', [45e3, 48e3, 73e3, 76e3], 'xticklabel', {'f_{p1}', 'f_{s1}', 'f_{s2}', 'f_{p2}'});
set(gca, 'YTick', [0.15, 0.85, 1, 1.15], 'yticklabel', {'\delta_2 = 0.15', '1 - \delta_1 = 0.85', '1', '1 + \delta_1 = 1.15'});
fvtool(Nz, Dz, 'Analysis', 'Phase');