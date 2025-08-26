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
objWaitBar = waitbar(0, "calcul en cours - "+num2str(0,1)+"%");
number_of_steps = 19;
current_step = 0;

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

syms m_face g I_xx_face I_yy_face I_zz_face
syms x(t) y(t) z(t) roll(t) pitch(t) yaw(t) sympi
syms alpha_arm_var(t) beta_arm k_tc
syms a_face d_face h_face
syms offset_prop

alpha_face = a_face/2 - d_face;
l0_face = (a_face + h_face)/2;
beta_face = (2*l0_face^2+6*alpha_face^2)^(1/2);
l_face = alpha_face*sin(alpha_arm_var)+(1/2*beta_face^2+alpha_face^2*(sin(alpha_arm_var))^2-alpha_face^2*(1+2*(cos(alpha_arm_var))^2))^(1/2);
m_drone = 6*m_face;

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

rot_BFi{1} = [ z0, -y0,  x0];
rot_BFi{2} = [ z0,  x0,  y0];
rot_BFi{3} = [ z0,  y0, -x0];
rot_BFi{4} = [ z0, -x0, -y0];
rot_BFi{5} = [ x0,  y0,  z0];
rot_BFi{6} = [ x0, -y0, -z0];
R_BF1 = rot_BFi{1}*Func_rot_z(alpha_arm_var);
R_BF2 = rot_BFi{2}*Func_rot_z(alpha_arm_var);
R_BF3 = rot_BFi{3}*Func_rot_z(alpha_arm_var);
R_BF4 = rot_BFi{4}*Func_rot_z(alpha_arm_var);
R_BF5 = rot_BFi{5}*Func_rot_z(alpha_arm_var);
R_BF6 = rot_BFi{6}*Func_rot_z(alpha_arm_var);

xb = R_0B*x0;
yb = R_0B*y0;
zb = R_0B*z0;

xF1 = R_0B*R_BF1*x0;
yF1 = R_0B*R_BF1*y0;
zF1 = R_0B*R_BF1*z0;
xF2 = R_0B*R_BF2*x0;
yF2 = R_0B*R_BF2*y0;
zF2 = R_0B*R_BF2*z0;
xF3 = R_0B*R_BF3*x0;
yF3 = R_0B*R_BF3*y0;
zF3 = R_0B*R_BF3*z0;
xF4 = R_0B*R_BF4*x0;
yF4 = R_0B*R_BF4*y0;
zF4 = R_0B*R_BF4*z0;
xF5 = R_0B*R_BF5*x0;
yF5 = R_0B*R_BF5*y0;
zF5 = R_0B*R_BF5*z0;
xF6 = R_0B*R_BF6*x0;
yF6 = R_0B*R_BF6*y0;
zF6 = R_0B*R_BF6*z0;

pF1 = [l_face; 0; 0];
pF2 = [0; l_face; 0];
pF3 = [-l_face;0; 0];
pF4 = [0;-l_face; 0];
pF5 = [0; 0; l_face];
pF6 = [0; 0;-l_face];

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

I_face = diag([I_xx_face + m_face*l_face^2, I_yy_face + m_face*l_face^2, I_zz_face]);

I_drone = R_BF1*I_face*R_BF1.' + R_BF2*I_face*R_BF2.' + R_BF3*I_face*R_BF3.' + R_BF4*I_face*R_BF4.' + R_BF5*I_face*R_BF5.' + R_BF6*I_face*R_BF6.';

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
T = 1/2*dq.'*M*dq;

U = -g_vec.'*m_drone*p;

L = T-U;

Tau_eq = diff(gradient(L, dq), t) - gradient(L, q);

Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

%% Identifications

text = 'Identification des matrices';
Func_Progress_Log(false, text, LOG_SIZE);

Mt = jacobian(Tau_eq, ddq);

C = Tau_eq - Mt*ddq;

Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

%% Passage de fonctions symboliques à équations symboliques

text = 'Substitutions des variables';
Func_Progress_Log(false, text, LOG_SIZE);


syms F1 F2 F3 F4 F5 F6 y1 y2 y3 y4 y5 y6 dy1 dy2 dy3 dy4 dy5 dy6 ddy1 ddy2 ddy3 ddy4 ddy5 ddy6 alpha_arm_val

y = [y1; y2; y3; y4; y5; y6];
dy = [dy1; dy2; dy3; dy4; dy5; dy6];
ddy = [ddy1; ddy2; ddy3; ddy4; ddy5; ddy6];

Tau_eq  =  subs(Tau_eq, [ddq; dq; q; alpha_arm_var], [ddy; dy; y; alpha_arm_val]);
Mt      =      subs(Mt, [ddq; dq; q; alpha_arm_var], [ddy; dy; y; alpha_arm_val]);
C       =       subs(C, [ddq; dq; q; alpha_arm_var], [ddy; dy; y; alpha_arm_val]);
S_eta   =   subs(S_eta, [ddq; dq; q; alpha_arm_var], [ddy; dy; y; alpha_arm_val]);
dS_eta  =  subs(dS_eta, [ddq; dq; q; alpha_arm_var], [ddy; dy; y; alpha_arm_val]);
l_face  =  subs(l_face, [ddq; dq; q; alpha_arm_var], [ddy; dy; y; alpha_arm_val]);
R_0B    =    subs(R_0B, [ddq; dq; q; alpha_arm_var], [ddy; dy; y; alpha_arm_val]);
R_BF1   =   subs(R_BF1, [ddq; dq; q; alpha_arm_var], [ddy; dy; y; alpha_arm_val]);
R_BF2   =   subs(R_BF2, [ddq; dq; q; alpha_arm_var], [ddy; dy; y; alpha_arm_val]);
R_BF3   =   subs(R_BF3, [ddq; dq; q; alpha_arm_var], [ddy; dy; y; alpha_arm_val]);
R_BF4   =   subs(R_BF4, [ddq; dq; q; alpha_arm_var], [ddy; dy; y; alpha_arm_val]);
pF1     =     subs(pF1, [ddq; dq; q; alpha_arm_var], [ddy; dy; y; alpha_arm_val]);
pF2     =     subs(pF2, [ddq; dq; q; alpha_arm_var], [ddy; dy; y; alpha_arm_val]);
pF3     =     subs(pF3, [ddq; dq; q; alpha_arm_var], [ddy; dy; y; alpha_arm_val]);
pF4     =     subs(pF4, [ddq; dq; q; alpha_arm_var], [ddy; dy; y; alpha_arm_val]);
I_drone = subs(I_drone, [ddq; dq; q; alpha_arm_var], [ddy; dy; y; alpha_arm_val]);

Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

%% Architecture du drone
%               ^
%               | F1
%               |
%               |
%    F2 <-------|-------> F4
%               |
%               |
%               | F3
%               -
%
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

p_rotor = sym('p_rotor_%d%d', [3,8]);
R_rotor = sym('R_rotor_%d%d%d', [3,3,8]);
R_rotor_viz = sym('R_rotor_viz_%d%d%d', [3,3,8]);
J_rotor_full = sym('J_rotor_%d%d', [6,8]);
J_rotor_flat_closed = sym('J_rotor_%d%d', [6,4]);
J_rotor_flat_opened = sym('J_rotor_%d%d', [6,4]);

p_rotor(:,1) = pF1 + R_BF1*[offset_prop; 0; h_face/2];
p_rotor(:,2) = pF2 + R_BF2*[offset_prop; 0; h_face/2];
p_rotor(:,3) = pF3 + R_BF3*[offset_prop; 0; h_face/2];
p_rotor(:,4) = pF4 + R_BF4*[offset_prop; 0; h_face/2];

p_rotor(:,5) = pF1 + R_BF1*[0; -offset_prop; h_face/2];
p_rotor(:,6) = pF2 + R_BF2*[0; -offset_prop; h_face/2];
p_rotor(:,7) = pF3 + R_BF3*[0; -offset_prop; h_face/2];
p_rotor(:,8) = pF4 + R_BF4*[0; -offset_prop; h_face/2];


R_rotor(:,:,1) = R_BF1*[-z0, y0, x0];
R_rotor(:,:,2) = R_BF2*[-z0, y0, x0];
R_rotor(:,:,3) = R_BF3*[-z0, y0, x0];
R_rotor(:,:,4) = R_BF4*[-z0, y0, x0];
R_rotor(:,:,5) = R_BF1*[x0, z0, -y0];
R_rotor(:,:,6) = R_BF2*[x0, z0, -y0];
R_rotor(:,:,7) = R_BF3*[x0, z0, -y0];
R_rotor(:,:,8) = R_BF4*[x0, z0, -y0];
 
for i = 1:8
    J_rotor_full(1:3,i) = R_rotor(:,:,i)*[0; 0; 1];
    if i<=4
        J_rotor_full(4:6,i) = R_rotor(:,:,i)*[0; 0; (-1)^(i)*k_tc] + cross(p_rotor(:,i), R_rotor(:,:,i)*[0; 0; 1]);
    else
        J_rotor_full(4:6,i) = R_rotor(:,:,i)*[0; 0; (-1)^(i-1)*k_tc] + cross(p_rotor(:,i), R_rotor(:,:,i)*[0; 0; 1]);
    end
end
for i = 1:4
    J_rotor_flat_closed(1:3,i) = [0; 0; 1];
    J_rotor_flat_closed(4:6,i) = [0; 0; (-1)^(i)*k_tc] + cross(p_rotor(:,i), [0; 0; 1]);
    J_rotor_flat_opened(1:3,i) = [0; 0; 1];
    J_rotor_flat_opened(4:6,i) = [0; 0; (-1)^(i+4)*k_tc] + cross(p_rotor(:,i+4), [0; 0; 1]);
end


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
m_face_num = m_drone_num/6;
a_face_num = 0.3;
d_face_num = 0.05;
h_face_num = 0.01;
offset_prop_num = a_face_num/4;
alpha_face_num = double(subs(alpha_face, [a_face, d_face], [a_face_num, d_face_num]));
beta_face_num = double(subs(beta_face, [a_face, d_face, h_face], [a_face_num, d_face_num, h_face_num]));


% alpha_arm_num = 30*pi/180;
alpha_arm_min = 0*pi/180;
alpha_arm_max = 90*pi/180;
alpha_arm_moy = (alpha_arm_min+alpha_arm_max)/2;
beta_arm_num = 0*pi/180;
dead_zone_alpha = 60*pi/180;

k_t_num = 0.1;
k_tc_num = 0.11;

fprop_sat_min = 0;
fprop_sat_max = 2/4*m_drone_num*g_ref;
prop_sat_min = fprop_sat_min/k_t_num;
prop_sat_max = fprop_sat_max/k_t_num;
fprop_sat_vec = [fprop_sat_min, fprop_sat_max];
prop_sat_vec = [prop_sat_min, prop_sat_max];

I_xx_num = 1/12*m_face_num*(3*a_face_num^2+h_face_num^2);
I_yy_num = 1/12*m_face_num*(3*a_face_num^2+h_face_num^2);
I_zz_num = 1/2*m_face_num*a_face_num^2;

% Remplacement dans les équations symboliques
replaced_syms = [g; sympi; m_face; a_face; h_face; offset_prop; d_face; beta_arm; k_tc; I_xx_face; I_yy_face; I_zz_face];
replacing_syms = [g_num; pi_num; m_face_num; a_face_num; h_face_num; offset_prop_num; d_face_num; beta_arm_num; k_tc_num; I_xx_num; I_yy_num; I_zz_num];

% Sym en fonction de y, alpha
S_eta_num = simp(formula(subs(S_eta, replaced_syms, replacing_syms)));
dS_eta_num = simp(formula(subs(dS_eta, replaced_syms, replacing_syms)));
I_drone_num = simp(formula(subs(I_drone, replaced_syms, replacing_syms)));
Mt_num = simp(formula(subs(Mt, replaced_syms, replacing_syms)));
C_num = simp(formula(subs(C, replaced_syms, replacing_syms)));
Tau_eq_num = simp(formula(subs(Tau_eq, replaced_syms, replacing_syms)));
l_face_num = simp(formula(subs(l_face, replaced_syms, replacing_syms)));
p_rotor_num = simp(formula(subs(p_rotor, replaced_syms, replacing_syms)));
J_rotor_full_num = simp(formula(subs(J_rotor_full, replaced_syms, replacing_syms)));
J_rotor_flat_closed_num = simp(formula(subs(J_rotor_flat_closed, replaced_syms, replacing_syms)));
J_rotor_flat_opened_num = simp(formula(subs(J_rotor_flat_opened, replaced_syms, replacing_syms)));
R_rotor_num = simp(formula(subs(R_rotor, replaced_syms, replacing_syms)));

% Variables purement numériques
% /

Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

%% Génération des fonctions

text = 'Génération de la fonction n°1';
Func_Progress_Log(false, text, LOG_SIZE);
matlabFunction(Mt_num, C_num, 'File', './Auto_Functions/AutoFunc_sym_function_drone', 'Vars', {y, dy, alpha_arm_val});
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
matlabFunction(J_rotor_full_num, J_rotor_flat_closed_num, J_rotor_flat_opened_num, 'File', './Auto_Functions/AutoFunc_sym_Mixer', 'Vars', {alpha_arm_val});
Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

text = 'Génération de la fonction n°4';
Func_Progress_Log(false, text, LOG_SIZE);
matlabFunction(l_face_num, 'File', './Auto_Functions/AutoFunc_sym_l_face', 'Vars', {alpha_arm_val});
Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

text = 'Génération de la fonction n°5';
Func_Progress_Log(false, text, LOG_SIZE);
matlabFunction(I_drone_num, 'File', './Auto_Functions/AutoFunc_sym_I_drone', 'Vars', {alpha_arm_val});
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

alpha0_num = 0;
l_face0_num = AutoFunc_sym_l_face(alpha0_num);

Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

%% Définition des visuels

text = 'Définition des visuels';
Func_Progress_Log(false, text, LOG_SIZE);

size_arrow = 20;
struct_arrow = size_arrow*[2, 0; 0.5, 1; 0.5, 0.5; -1.5, 0.5; -1.5, -0.5; 0.5, -0.5; 0.5, -1;];
h_arrow = size_arrow*0.5;

face_color = [0.0, 0.8, 0.8;
              0.8, 0.0, 0.5;
              0.8, 0.8, 0.0;
              0.0, 0.8, 0.0;
              0.8, 0.8, 0.8;
              0.2, 0.2, 0.2];
face_opacity = 1;
triangle_color = [0.5, 0.5, 0.5];
triangle_opacity = 1;
rotor_color_x = [1 0.4118 0.1608];
rotor_color_y = [0.4118 0.1608 1];
rotor_opacity = 1;

Lims = [-10, 10; -10, 10; -2, 35; -20, 20];

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

Script_PreCompute_SmartCascade_Juno;

%% Définition de la trajectoire

text = 'Définition de la trajectoire';
Func_Progress_Log(false, text, LOG_SIZE);

StopTime = 5;

% Time points
Trajectory_time = [0 2.5 19 21 37.5 40];
Attitude_time = Trajectory_time;
Amplitude_sinus = 0*1.5;
Pulsation_sinus = 2*pi/4;

% p(t)
Trajectory_points = [   init(1) 0.5 0.5 0.5 0.5 0
                        init(2) 0 0 1 1 1
                        init(3) 1 1 1 1 1];

% 1dp(t)
Trajectory_velBC = [    init(7) 0 0 0 0 0
                        init(8) 0 0 0 0 0
                        init(9) 0 0 0 0 0];

% 2dp(t)
Trajectory_accelBC = zeros(size(Trajectory_points));


% 3dp(t)
Trajectory_jerkBC = zeros(size(Trajectory_points));

% 4dp(t)
Trajectory_snapBC = zeros(size(Trajectory_points));

% eta(t)
Attitude_points = [     init(4) 0 0 0 0 0
                        init(5) 0 0 0 0 0
                        init(6) 0 0 0 0 0];

% 1d_eta(t)
Attitude_velBC = [      init(7) 0 0 0 0 0
                        init(8) 0 0 0 0 0
                        init(9) 0 0 0 0 0];

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