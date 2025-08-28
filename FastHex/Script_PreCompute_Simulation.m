%% INIT

clear variables
close all
addpath("Functions\")
addpath("Auto_Functions\")
addpath("Images\")
addpath("Videos\")
fprintf('_____________________________________________________________________\n')
fprintf('Nouveau cacul........................................................\n')
fprintf('_____________________________________________________________________\n')
number_of_steps = 17;
current_step = 0;
objWaitBar = waitbar(current_step/number_of_steps, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%");

%% DEFINITIONS

LOG_SIZE = 60;
text = 'Définitions fonctions locales';
Func_Progress_Log(false, text, LOG_SIZE);

% Fonction pour sélectionner un élément dans un vecteur (taille 3) de fonctions sans évaluer la fonction
id3 = @(vec, i) vec.'*[abs(abs(sign(i-1))-1); abs(abs(sign(i-2))-1); abs(abs(sign(i-3))-1)];

% Opérateur chapeau
skew = @(x) [0, -id3(x,3), id3(x,2); id3(x,3), 0, -id3(x,1); -id3(x,2), id3(x,1), 0];

% Création d'une matrice de rotation selon z0 pour R²
rot_z2 = @(q) [cos(q) -sin(q)
               sin(q)  cos(q)];

% Poignée plus courte pour le simplify()
simp = @(x) simplify(x);

% Calcul de coefficient de détermination
R2 = @(ts1, ts2) 1 - sum((ts1 - ts2).^2) / sum((ts1 - mean(ts2)).^2);

% Axes du repère initial
x0 = [1;0;0];
y0 = [0;1;0];
z0 = [0;0;1];

Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

%% Variables symboliques

text = 'Défintion variables symboliques';
Func_Progress_Log(false, text, LOG_SIZE);

syms m_drone g I_xx I_yy I_zz
syms x(t) y(t) z(t) roll(t) pitch(t) yaw(t) sympi 
syms h_drone r_arm alpha_arm beta_arm k_tc

% Dérivées des variables d'état
dx  = diff(x, t, 1);
dy  = diff(y, t, 1);
dz  = diff(z, t, 1);
ddx = diff(x, t, 2);
ddy = diff(y, t, 2);
ddz = diff(z, t, 2);
p = [x; y; z];
dp = [dx; dy; dz];
ddp = [ddx; ddy; ddz];

eta    = [roll; pitch; yaw];
d_eta  = diff(eta, t, 1);
dd_eta = diff(eta, t, 2);

R_0B = Func_rot_z(yaw)*Func_rot_y(pitch)*Func_rot_x(roll);

xb = R_0B*[1; 0; 0];
yb = R_0B*[0; 1; 0];
zb = R_0B*[0; 0; 1];

S_eta = [1       0         -sin(pitch)    ;
         0  cos(roll) sin(roll)*cos(pitch);
         0 -sin(roll) cos(roll)*cos(pitch)];

dS_eta = diff(S_eta, t);

omega   = S_eta*d_eta;
domega  = diff(omega, t);

% Et vecteurs associés
q = [p; eta];
dq = [dp; d_eta];
ddq = [ddp; dd_eta];

% Matrice d'inertie du drone

I_drone = diag([I_xx, I_yy, I_zz]);

% Vecteur de gravité
g_vec = [0; 0; g];

Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

%% Calculs intermédiaires

text = 'Calculs intermédiaires';
Func_Progress_Log(false, text, LOG_SIZE);

Mb = [m_drone*eye(3), zeros(3); zeros(3), I_drone];

Jt = [R_0B.', zeros(3); zeros(3), S_eta];

M = Jt.'*Mb*Jt;

Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

%% Calcul des Lagrangiens

text = 'Calculs des Lagrangiens';
Func_Progress_Log(false, text, LOG_SIZE);

% T = 1/2*m_drone*(dp.'*dp) + 1/2*omega.'*I_drone*omega;
T = simp(1/2*dq.'*M*dq);

U = simp(-g_vec.'*m_drone*p);

L = T-U;

Tau_eq = simp(diff(gradient(L, dq), t) - gradient(L, q));

Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

%% Identifications

text = 'Identification des matrices';
Func_Progress_Log(false, text, LOG_SIZE);

Mt = simp(jacobian(Tau_eq, ddq));

C = simp(Tau_eq - Mt*ddq);

Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

%% Passage de fonctions symboliques à équations symboliques

text = 'Substitutions des variables';
Func_Progress_Log(false, text, LOG_SIZE);

syms F1 F2 F3 F4 F5 F6 y1 y2 y3 y4 y5 y6 dy1 dy2 dy3 dy4 dy5 dy6 ddy1 ddy2 ddy3 ddy4 ddy5 ddy6

y = [y1; y2; y3; y4; y5; y6];
dy = [dy1; dy2; dy3; dy4; dy5; dy6];
ddy = [ddy1; ddy2; ddy3; ddy4; ddy5; ddy6];

Tau_eq = subs(Tau_eq, [ddq; dq; q], [ddy; dy; y]);
Mt     =     subs(Mt, [ddq; dq; q], [ddy; dy; y]);
C      =      subs(C, [ddq; dq; q], [ddy; dy; y]);
R_0B   =   subs(R_0B, [ddq; dq; q], [ddy; dy; y]);
S_eta  =  subs(S_eta, [ddq; dq; q], [ddy; dy; y]);
dS_eta = subs(dS_eta, [ddq; dq; q], [ddy; dy; y]);

Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")


%% Architecture du drone
%    
%    6-/      \+5
%       \    /
%        \  /
%         \/
% 1+|------------|-4
%         /\
%        /  \
%       /    \
%    2-\      /+3
%
%               ^
%               | x_b
%               |
%               |
%   y_b <-------|
%   

text = ['Définition de l' char("'") 'architecture du drone'];
Func_Progress_Log(false, text, LOG_SIZE);

% Position et orientation des hélices

p_rotor = sym('p_rotor_%d%d', [3,6]);
R_rotor = sym('R_rotor_%d%d%d', [3,3,6]);
R_rotor_viz = sym('R_rotor_viz_%d%d%d', [3,3,6]);
J_rotor = sym('J_rotor_%d%d', [6,6]);
J_rotor_flat = sym('J_rotor_flat_%d%d', [6,6]);

alpha_struct = [1 -1 1 -1 1 -1];
beta_struct = [1 1 1 1 1 1];

for i = 1:6
    p_rotor(:,i) = r_arm*[cospi((i-1)/3+1/2); sinpi((i-1)/3+1/2); 0];
    R_rotor(:,:,i) = Func_rot_z((i-1)*sympi/3+sympi/2)*Func_rot_y(beta_struct(i)*beta_arm)*Func_rot_x(alpha_struct(i)*alpha_arm);
    R_rotor_viz(:,:,i) = Func_rot_z((i-1)*sympi/3+sympi/2)*Func_rot_y(beta_struct(i)*beta_arm);
end

for i = 1:6
    J_rotor(1:3,i) = R_rotor(:,:,i)*[0; 0; 1];
    J_rotor_flat(1:3,i) = [0; 0; 1];
    J_rotor(4:6,i) = R_rotor(:,:,i)*[0; 0; (-1)^(i)*k_tc] + cross(p_rotor(:,i), R_rotor(:,:,i)*[0; 0; 1]);
    J_rotor_flat(4:6,i) = [0; 0; (-1)^(i)*k_tc] + cross(p_rotor(:,i), [0; 0; 1]);
end

F_vec = [F1; F2; F3; F4; F5; F6];
Wrench = J_rotor*F_vec;

Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

%% Calcul numérique

text = 'Préparation des calculs numériques';
Func_Progress_Log(false, text, LOG_SIZE);

% Définition du temps
time_start = 0;
time_end = 30;
time_step = 1e-1;
time = time_start:time_step:time_end;

% Définition de la gravité et de pi
g_ref = 9.80665;
g_num = -g_ref;
pi_num = pi;
    
% Définition de constantes du drone

m_drone_num = 1;
r_arm_num = 0.35;
h_drone_num = 0.05;

% alpha_arm_num = 30*pi/180;
alpha_arm_min = 0*pi/180;
alpha_arm_max = 50*pi/180;
alpha_arm_moy = (alpha_arm_min+alpha_arm_max)/2;
beta_arm_num = 0*pi/180;
dead_zone_alpha = 1*pi/180;

k_t_num = 0.01;
k_tc_num = 0.11;

fprop_sat_min = 0;
fprop_sat_max = 2/6*m_drone_num*g_ref/cos(alpha_arm_moy);
prop_sat_min = fprop_sat_min/k_t_num;
prop_sat_max = fprop_sat_max/k_t_num;
fprop_sat_vec = [fprop_sat_min, fprop_sat_max];
prop_sat_vec = [prop_sat_min, prop_sat_max];
fprop_sat_vec_margin = [fprop_sat_min + 0.1*(fprop_sat_max - fprop_sat_min), fprop_sat_max - 0.1*(fprop_sat_max - fprop_sat_min)];

I_xx_num = 1/12*m_drone_num*(3*r_arm_num^2+h_drone_num^2);
I_yy_num = 1/12*m_drone_num*(3*r_arm_num^2+h_drone_num^2);
I_zz_num = 1/2*m_drone_num*r_arm_num^2;

% Remplacement dans les équations symboliques
% replaced_syms = [g; sympi; m_drone; r_arm; h_drone; alpha_arm; beta_arm; k_tc; I_xx; I_yy; I_zz];
% replacing_syms = [g_num; pi_num; m_drone_num; r_arm_num; h_drone_num; alpha_arm_num; beta_arm_num; k_tc_num; I_xx_num; I_yy_num; I_zz_num];
replaced_syms = [g; sympi; m_drone; r_arm; h_drone; beta_arm; k_tc; I_xx; I_yy; I_zz];
replacing_syms = [g_num; pi_num; m_drone_num; r_arm_num; h_drone_num; beta_arm_num; k_tc_num; I_xx_num; I_yy_num; I_zz_num];

% Sym en fonction de y, ...
Tau_eq_num = formula(subs(Tau_eq, replaced_syms, replacing_syms));
Wrench_num = formula(subs(Wrench, replaced_syms, replacing_syms));
S_eta_num = formula(subs(S_eta, replaced_syms, replacing_syms));
dS_eta_num = formula(subs(dS_eta, replaced_syms, replacing_syms));
Mt_num = formula(subs(Mt, replaced_syms, replacing_syms));
C_num = formula(subs(C, replaced_syms, replacing_syms));
J_rotor_num = formula(subs(J_rotor, replaced_syms, replacing_syms));
R_rotor_num = formula(subs(R_rotor, replaced_syms, replacing_syms));

% Variables purement numériques
I_drone_num = double(subs(I_drone, replaced_syms, replacing_syms));
p_rotor_num = double(subs(p_rotor, replaced_syms, replacing_syms));
R_rotor_viz_num = double(subs(R_rotor_viz, replaced_syms, replacing_syms));
J_rotor_flat_num = double(subs(J_rotor_flat, replaced_syms, replacing_syms));

I_moments = diag(I_drone_num);
I_products = zeros([3, 1]);
Jinv_rotor_num = pinv(J_rotor_num);

Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

%% Génération des fonctions

text = 'Génération de la fonction n°1';
Func_Progress_Log(false, text, LOG_SIZE);
matlabFunction(Mt_num, C_num, 'File', './Auto_Functions/AutoFunc_sym_function_drone', 'Vars', {y, dy, F_vec});
Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

text = 'Génération de la fonction n°2';
Func_Progress_Log(false, text, LOG_SIZE);
matlabFunction(S_eta_num, dS_eta_num, 'File', './Auto_Functions/AutoFunc_S_eta_gen', 'Vars', {y, dy});
Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

text = 'Génération de la fonction n°3';
Func_Progress_Log(false, text, LOG_SIZE);
matlabFunction(J_rotor_num, 'File', './Auto_Functions/AutoFunc_sym_Mixer', 'Vars', {alpha_arm});
Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

%% Définition des états initiaux

text = 'Défintion des états initiaux';
Func_Progress_Log(false, text, LOG_SIZE);

p0 = [0; 0; 1];
dp0 = [0; 0; 0];

eta0_XYZ = [0; 0; 0];
omega0 = [0; 0; 0];

R_init = Func_rot_z(eta0_XYZ(3))*Func_rot_y(eta0_XYZ(2))*Func_rot_x(eta0_XYZ(1));
eta0_xyz = rotm2eul(R_init, "XYZ").';

[S_eta0_XYZ] = AutoFunc_S_eta_gen([p0; eta0_XYZ], zeros([6,1]));
d_eta0_XYZ = S_eta0_XYZ\omega0;

[S_eta0_xyz] = AutoFunc_S_eta_gen([p0; eta0_xyz], zeros([6,1]));
d_eta0_xyz = S_eta0_xyz\omega0;

init         = [p0; eta0_XYZ; dp0; d_eta0_XYZ];
init_bushing = [p0; eta0_xyz; dp0; d_eta0_xyz];

Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

%% Définition des visuels

text = 'Définition des visuels';
Func_Progress_Log(false, text, LOG_SIZE);

size_arrow = 20;
struct_arrow = size_arrow*[2, 0; 0.5, 1; 0.5, 0.5; -1.5, 0.5; -1.5, -0.5; 0.5, -0.5; 0.5, -1;];
h_arrow = size_arrow*0.5;

h_hexa = h_drone_num;
L_hexa = r_arm_num;
r_hexa = h_hexa/sin(pi/3)+h_hexa*tan(pi/2-pi/3);
struct_bras_hexa = [h_hexa, r_hexa; h_hexa, L_hexa; -h_hexa, L_hexa; -h_hexa, r_hexa];

struct_hexa = zeros(19,2);
struct_hexa(1:4,:) = struct_bras_hexa;

for i = 1:5
    struct_hexa(4+i*3-2:4+i*3, :) = [(rot_z2(i*pi/3)*struct_bras_hexa(2,:).').'; (rot_z2(i*pi/3)*struct_bras_hexa(3,:).').'; (rot_z2(i*pi/3)*struct_bras_hexa(4,:).').'];
end

struct_hexa = struct_hexa(1:18, :);

Tile_delta = 0.5;
Tile_A = [211 204 204]/255;
Tile_B = [227 230 232]/255;

Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

%% Controller parameters

text = 'Paramètrage du controlleur';
Func_Progress_Log(false, text, LOG_SIZE);

k_corr_x = 5;
k_corr_v = 2.5;
k_corr_R = [16 0 0
    0 16 0
    0 0 30];

k_corr_Omega = [1.5 0  0
    0 1.5 0
    0  0 4];


Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

%% FastHex control verifications precompute

text = 'Calcul préliminaires du controlleur';
Func_Progress_Log(false, text, LOG_SIZE);

% Script_PreCompute_SmartCascade;

current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")
Func_Progress_Log(true, text, LOG_SIZE);

%% Définition de la trajectoire

text = 'Définition de la trajectoire';
Func_Progress_Log(false, text, LOG_SIZE);

% Time points
Trajectory_time = [0 2.5 5 7.5 10 12.5 15 17.5 20 22.5 25];
Attitude_time = Trajectory_time;
Amplitude_sinus = 1.4;
Pulsation_sinus = 2*pi/4;

% p0 = [-Amplitude_sinus; 0; 1];
% dp0 = [0; Amplitude_sinus*Pulsation_sinus; 0];
% init         = [p0; eta0_XYZ; dp0; d_eta0_XYZ];
% init_bushing = [p0; eta0_xyz; dp0; d_eta0_xyz];

% p(t)
Trajectory_points = [init(1) 2 4 0 -2 0 -4 0 0 0 0
    init(2) 0 1 0  2 0 -2 0 0 0 0
    init(3) 1 1 1  1 2  1 1 1 1 1];

% 1dp(t)
Trajectory_velBC = [init(7) 0 0 -6/5 0   0 0 0 0 0 0
    init(8) 0 0    0 0.5 0 0 0 0 0 0
    init(9) 0 0    0 0   0 0 0 0 0 0];

% 2dp(t)
Trajectory_accelBC = zeros(size(Trajectory_points));

% 3dp(t)
Trajectory_jerkBC = zeros(size(Trajectory_points));

% 4dp(t)
Trajectory_snapBC = zeros(size(Trajectory_points));

% eta(t)
Attitude_points = [init(4) pi/3 0    0 0    0 0 pi/3 0  0.4*pi 0
    init(5) pi/3 0 pi/4 0    0 0 pi/3 0 -0.4*pi 0
    init(6) pi/3 0    0 0 pi/4 0    0 0    pi/6 0];

% 1d_eta(t)
Attitude_velBC = [init(7) 0 0 0 0 0 0 0 0 0 0
    init(8) 0 0 0 0 0 0 0 0 0 0
    init(9) 0 0 0 0 0 0 0 0 0 0];

% 2d_eta(t)
Attitude_accelBC = zeros(size(Attitude_points));

% 3d_eta(t)
Attitude_jerkBC = zeros(size(Attitude_points));

% 4d_eta(t)
Attitude_snapBC = zeros(size(Attitude_points));

Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

%% FIN

fprintf('_____________________________________________________________________\n')
fprintf('Calcul prêt..........................................................\n')
fprintf('_____________________________________________________________________\n')
close(objWaitBar)