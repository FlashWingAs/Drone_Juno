%% Calculs symboliques

text = 'Pré-calculs symboliques du controlleur';
Func_Progress_Log(false, text, LOG_SIZE);

% Espace de forces moteurs pour un torseur donné
syms Ft Fz yaw alpha

Fx = cos(yaw)*Ft;
Fy = sin(yaw)*Ft;
Fb_test = [Fx; Fy; Fz; 0; 0; 0];
J_mot = AutoFunc_sym_Mixer(alpha);
F_mot = J_mot\Fb_test;
% F_mot_eq = formula(subs(F_mot, yaw, Yaw_setpoint));
F_mot_eq = formula(F_mot);

Fmot1 = F_mot_eq(1);
Fmot2 = F_mot_eq(2);
Fmot3 = F_mot_eq(3);

% Lien des efforts tangeantiels par les vitesses moteurs

syms w1 w2

Fz_w2 = solve(Fmot1==w1, Fz);
Fmot2_w2 = subs(Fmot2, Fz, Fz_w2);
Ft_w2 = solve(Fmot2_w2==w2, Ft);
Fz_w2 = subs(Fz_w2, Ft, Ft_w2);

Fz_w3 = solve(Fmot1==w1, Fz);
Fmot3_w3 = subs(Fmot3, Fz, Fz_w3);
Ft_w3 = solve(Fmot3_w3==w2, Ft);
Fz_w3 = subs(Fz_w3, Ft, Ft_w3);

F_mot_eq_w2 = simplify(subs(F_mot_eq, [Ft, Fz], [Ft_w2, Fz_w2]));
Fb_test_w2 = J_mot*F_mot_eq_w2;
k3x_w2 = diff(F_mot_eq_w2(3), w1);
k3y_w2 = diff(F_mot_eq_w2(3), w2);
c3_w2 = F_mot_eq_w2(3)-(k3x_w2*w1+k3y_w2*w2);


F_mot_eq_w3 = simplify(subs(F_mot_eq, [Ft, Fz], [Ft_w3, Fz_w3]));
Fb_test_w3 = J_mot*F_mot_eq_w3;
k3x_w3 = diff(F_mot_eq_w3(3), w1);
k3y_w3 = diff(F_mot_eq_w3(3), w2);
c3_w3 = F_mot_eq_w3(3)-(k3x_w3*w1+k3y_w3*w2);

Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

%% Génération des fonctions

text = 'Génération des fonctions internes du controlleur';
Func_Progress_Log(false, text, LOG_SIZE);

matlabFunction(Fb_test_w2, 'File', './Auto_Functions/AutoFunc_sym_Fb_trans_w2', 'Vars', {alpha, yaw, w1, w2});
matlabFunction(F_mot_eq_w2(3), 'File', './Auto_Functions/AutoFunc_sym_F3_trans_w2', 'Vars', {alpha, yaw, w1, w2});
matlabFunction(k3x_w2, 'File', './Auto_Functions/AutoFunc_sym_k3x_w2', 'Vars', {alpha, yaw, w1, w2});
matlabFunction(k3y_w2, 'File', './Auto_Functions/AutoFunc_sym_k3y_w2', 'Vars', {alpha, yaw, w1, w2});
matlabFunction(c3_w2, 'File', './Auto_Functions/AutoFunc_sym_c3_w2', 'Vars', {alpha, yaw, w1, w2});

matlabFunction(Fb_test_w3, 'File', './Auto_Functions/AutoFunc_sym_Fb_trans_w3', 'Vars', {alpha, yaw, w1, w2});
matlabFunction(F_mot_eq_w3(3), 'File', './Auto_Functions/AutoFunc_sym_F3_trans_w3', 'Vars', {alpha, yaw, w1, w2});
matlabFunction(k3x_w3, 'File', './Auto_Functions/AutoFunc_sym_k3x_w3', 'Vars', {alpha, yaw, w1, w2});
matlabFunction(k3y_w3, 'File', './Auto_Functions/AutoFunc_sym_k3y_w3', 'Vars', {alpha, yaw, w1, w2});
matlabFunction(c3_w3, 'File', './Auto_Functions/AutoFunc_sym_c3_w3', 'Vars', {alpha, yaw, w1, w2});

Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")