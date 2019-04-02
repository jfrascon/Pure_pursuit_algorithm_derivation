clc
close all;
clear all;

rb = 0.35;
n = 1.5;
N = 25;
yaw_rot = (60.0/180.0)*pi;
yaw_rotmax = 0.5*pi;
Vmax = 10*1000/3600;
Wmax = (60/180)*pi;
Wmaxrot = (60/180)*pi;
Wminrot = Wmax/2;

nsamples = 100;
nsamples_bot = nsamples/2 - 15;
nsamples_top = nsamples/2 + 15;


lah = n*(2*rb);
yaw_min = (10.0/180.0)*pi;

% -------------------------------



yaw_pi2rotmax = linspace(-pi, -yaw_rotmax, nsamples);
yaw_rotmax2rot = linspace(-yaw_rotmax, -yaw_rot, nsamples);
yaw_rot2min = linspace(-yaw_rot, -yaw_min, nsamples);
yaw_min2zer = linspace(-yaw_min, 0, nsamples);
yaw_zer2min = linspace(0, yaw_min, nsamples);
yaw_min2rot = linspace(yaw_min, yaw_rot, nsamples);
yaw_rot2rotmax = linspace(yaw_rot, yaw_rotmax, nsamples);
yaw_rotmax2pi = linspace(yaw_rotmax, pi, nsamples);

yaw = [yaw_pi2rotmax yaw_rotmax2rot yaw_rot2min yaw_min2zer yaw_zer2min yaw_min2rot yaw_rot2rotmax yaw_rotmax2pi];
yaw_deg = (180.0/pi)*yaw;


C_pi2rotmax = 2*sin(yaw_pi2rotmax)/lah;
C_rotmax2rot = 2*sin(yaw_rotmax2rot)/lah;
C_rot2min = 2*sin(yaw_rot2min)/lah;
C_min2zer = 2*sin(yaw_min2zer)/lah;
C_zer2min = 2*sin(yaw_zer2min)/lah;
C_min2rot = 2*sin(yaw_min2rot)/lah;
C_rot2rotmax = 2*sin(yaw_rot2rotmax)/lah;
C_rotmax2pi = 2*sin(yaw_rotmax2pi)/lah;

C = [C_rot2min C_min2zer C_zer2min C_min2rot];

figure(1)
plot(180/pi*[yaw_rot2min yaw_min2zer yaw_zer2min yaw_min2rot], C, 'c', 'linewidth', 2);
grid on;

V_pi2rotmax = zeros(size(yaw_pi2rotmax));
V_rotmax2rot = zeros(size(yaw_rotmax2rot));
V_rot2min = (Vmax/(yaw_rot-yaw_min)).*(yaw_rot2min+yaw_min)+Vmax;
V_min2zer = Vmax*ones(size(yaw_min2zer));
V_zer2min = Vmax*ones(size(yaw_zer2min));
V_min2rot = (-Vmax/(yaw_rot-yaw_min)).*(yaw_min2rot-yaw_min)+Vmax;
V_rot2rotmax = zeros(size(yaw_rot2rotmax));
V_rotmax2pi = zeros(size(yaw_rotmax2pi));

V = [V_pi2rotmax V_rotmax2rot V_rot2min V_min2zer V_zer2min V_min2rot V_rot2rotmax V_rotmax2pi];

figure(2);
plot(yaw_deg, V,  'r', 'linewidth', 2);
yyaxis left
grid on;

W_pi2rotmax = -Wmaxrot*ones(size(yaw_pi2rotmax));
W_rotmax2rot = (Wmaxrot-Wminrot)/(yaw_rotmax-yaw_rot)*(yaw_rotmax2rot+yaw_rot)-Wminrot;
W_rot2min = V_rot2min.*C_rot2min;
W_min2zer = zeros(size(V_min2zer));
W_zer2min = zeros(size(V_min2zer));
W_min2rot = V_min2rot.*C_min2rot;
W_rot2rotmax = (Wmaxrot-Wminrot)/(yaw_rotmax-yaw_rot)*(yaw_rot2rotmax-yaw_rot)+Wminrot;
W_rotmax2pi = Wmaxrot*ones(size(yaw_rotmax2pi));

W = [W_pi2rotmax W_rotmax2rot W_rot2min W_min2zer W_zer2min W_min2rot W_rot2rotmax W_rotmax2pi];

figure(2);
yyaxis right
plot(yaw_deg, W, 'g', 'linewidth', 2);
hold on
grid on;

W_rot2min(W_rot2min < -Wmax)= -Wmax;
W_min2rot(W_min2rot>Wmax) = Wmax;

W_trunc = [W_pi2rotmax W_rotmax2rot W_rot2min W_min2zer W_zer2min W_min2rot W_rot2rotmax W_rotmax2pi];

figure(2)
yyaxis right
plot(yaw_deg, W_trunc, 'm.-', 'linewidth', 2);

V_rot2min = W_rot2min./C_rot2min;
V_min2rot = W_min2rot./C_min2rot;

V_trunc = [V_pi2rotmax V_rotmax2rot V_rot2min V_min2zer V_zer2min V_min2rot V_rot2rotmax V_rotmax2pi];

figure(2)
hold on;
yyaxis left
plot(yaw_deg, V_trunc, 'b.', 'linewidth', 2);

figure(3)
yyaxis left
plot(180/pi*[yaw_rot2min yaw_min2zer yaw_zer2min yaw_min2rot], C, 'c', 'linewidth', 2);
hold on;
yyaxis right
plot(yaw_deg, V,  'r', 'linewidth', 2);
yyaxis right
plot(yaw_deg, V_trunc, 'b.-', 'linewidth', 2);
grid on;

figure(4)
yyaxis left
plot(180/pi*[yaw_rot2min yaw_min2zer yaw_zer2min yaw_min2rot], C, 'c', 'linewidth', 2);
hold on;
yyaxis right
plot(yaw_deg, W,  'g', 'linewidth', 2);
yyaxis right
plot(yaw_deg, W_trunc, 'm.-', 'linewidth', 2);
grid on;

%{
yaw_l3 = linspace(-0.95*pi, -yaw_rot_max, nsamples);
yaw_l2 = linspace(-yaw_rot_max, -yaw_rot, nsamples);
yaw_l1 = linspace(-yaw_rot, -yaw_min, nsamples);

yaw_c = linspace(-yaw_min, yaw_min, nsamples);

yaw_r1 = linspace(yaw_min, yaw_rot, nsamples);
yaw_r2 = linspace(yaw_rot, yaw_rot_max, nsamples);
yaw_r3 = linspace(yaw_rot_max, 0.95*pi, nsamples);

yaw_l3_deg = yaw_l3 * 180/pi;
yaw_l2_deg = yaw_l2 * 180/pi;
yaw_l1_deg = yaw_l1 * 180/pi;

yaw_c_deg = yaw_c * 180/pi;

yaw_r1_deg = yaw_r1 * 180/pi;
yaw_r2_deg = yaw_r2 * 180/pi;
yaw_r3_deg = yaw_r3 * 180/pi;

Cs_l3_theor = lah./(2*sin(yaw_l3));
Cs_l2_theor = lah./(2*sin(yaw_l2));
s_l1_theor = lah./(2*sin(yaw_l1));

Rs_lc_theor = lah./(2*sin(yaw_c(1:nsamples_bot)));
Rs_rc_theor = lah./(2*sin(yaw_c(nsamples_top:end)));

Rs_r1_theor = lah./(2*sin(yaw_r1));
Rs_r2_theor = lah./(2*sin(yaw_r2));
Rs_r3_theor = lah./(2*sin(yaw_r3));

figure(1)
plot(yaw_l3_deg, Rs_l3_theor, 'r-');
hold on
plot(yaw_l2_deg, Rs_l2_theor, 'r-');
plot(yaw_l1_deg, Rs_l1_theor, 'r-');

plot(yaw_c_deg(1:nsamples_bot), Rs_lc_theor, 'g-');
plot(yaw_c_deg(nsamples_top:end), Rs_rc_theor, 'g-');

plot(yaw_r1_deg, Rs_r1_theor, 'b-');
plot(yaw_r2_deg, Rs_r2_theor, 'b-');
plot(yaw_r3_deg, Rs_r3_theor, 'b-');
grid()

Rs_l3 = zeros(size(yaw_l3));
Rs_l2 = zeros(size(yaw_l2));
Rs_l1 = lah./(2*sin(yaw_l1));

%Rs_lc = -Rmax*ones(size(yaw_lc));
%Rs_rc = Rmax*ones(size(yaw_rc));

Rs_r1 = lah./(2*sin(yaw_r1));
Rs_r2 = zeros(size(yaw_r2));
Rs_r3 = zeros(size(yaw_r3));

plot(yaw_l3_deg, Rs_l3, 'c', 'linewidth', 2);
plot(yaw_l2_deg, Rs_l2, 'c', 'linewidth', 2);
plot(yaw_l1_deg, Rs_l1, 'c', 'linewidth', 2);

plot(yaw_r1_deg, Rs_r1, 'c', 'linewidth', 2);
plot(yaw_r2_deg, Rs_r2, 'c', 'linewidth', 2);
plot(yaw_r3_deg, Rs_r3, 'c', 'linewidth', 2);

% -------------------------------

V_l3 = zeros(size(yaw_l3));
V_l2 = zeros(size(yaw_l2));
m = -Vmax/(yaw_rot-yaw_min);
V_l1 = m*(abs(yaw_l1)-yaw_min) + Vmax;

V_c = Vmax*ones(size(yaw_c));

V_r1 = m*(abs(yaw_r1)-yaw_min) + Vmax;
V_r2 = zeros(size(yaw_r2));
V_r3 = zeros(size(yaw_r3));

figure(2)
plot(yaw_l3_deg, V_l3, 'c', 'linewidth', 2);
hold on;
grid on;
plot(yaw_l2_deg, V_l2, 'c', 'linewidth', 2);
plot(yaw_l1_deg, V_l1, 'c', 'linewidth', 2);

plot(yaw_c_deg, V_c, 'c', 'linewidth', 2);

plot(yaw_r1_deg, V_r1, 'c', 'linewidth', 2);
plot(yaw_r2_deg, V_r2, 'c', 'linewidth', 2);
plot(yaw_r3_deg, V_r3, 'c', 'linewidth', 2);

% -------------------------------

mrot = (Wmaxrot-Wminrot)/(yaw_rot_max-yaw_rot);

W_l3 = -Wmaxrot*ones(size(yaw_l3));
W_l2 = mrot*(yaw_l2 + yaw_rot) - Wminrot;
W_l1 = V_l1./Rs_l1;

W_c = zeros(size(yaw_c));

W_r1 = V_r1./Rs_r1;
W_r2 = mrot*(yaw_r2 - yaw_rot) + Wminrot;
W_r3 = Wmaxrot*ones(size(yaw_r3));
yaw_deg
figure(3)
plot(yaw_l3_deg, W_l3, 'c', 'linewidth', 2);
hold on;
grid on;
plot(yaw_l2_deg, W_l2, 'c', 'linewidth', 2);
plot(yaw_l1_deg, W_l1, 'c', 'linewidth', 2);

plot(yaw_c_deg, W_c, 'c', 'linewidth', 2);

plot(yaw_r1_deg, W_r1, 'c', 'linewidth', 2);
plot(yaw_r2_deg, W_r2, 'c', 'linewidth', 2);
plot(yaw_r3_deg, W_r3, 'c', 'linewidth', 2);

W_l1(W_l1 < -Wmax) = -Wmax;

W_r1(W_r1 > Wmax) = Wmax;

plot(yaw_l1_deg, W_l1, 'b--', 'linewidth', 2);
plot(yaw_r1_deg, W_r1, 'b--', 'linewidth', 2);

V_l1 = W_l1.* Rs_l1; 
V_r1 = W_r1.* Rs_r1;

figure(2);
plot(yaw_l1_deg, V_l1, 'b--', 'linewidth', 2);
plot(yaw_r1_deg, V_r1, 'b--', 'linewidth', 2);
%}

