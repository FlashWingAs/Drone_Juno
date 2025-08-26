% Espace de forces moteurs pour un torseur donné
syms Fx Fy Fz alpha

Fb_test = [Fx; Fy; Fz; 0; 0; 0];
J_mot = AutoFunc_sym_Mixer(alpha);
F_mot = J_mot\Fb_test;
% F_mot_eq = formula(subs(F_mot, yaw, Yaw_setpoint));
F_mot_eq = formula(F_mot);

Fmot1 = F_mot_eq(1);
Fmot2 = F_mot_eq(2);

% Lien des efforts tangeantiels par les vitesses moteurs

syms w1 w2

Fx_w = solve(Fmot1==w1, Fx);
Fmot2_w = subs(Fmot2, Fx, Fx_w);
Fy_w = solve(Fmot2_w==w2, Fy);

F_mot_eq_w = simplify(subs(F_mot_eq, [Fx, Fy], [Fx_w, Fy_w]));
Fb_test_w = J_mot*F_mot_eq_w;
k3x = diff(F_mot_eq_w(3), w1);
k3y = diff(F_mot_eq_w(3), w2);
c3 = F_mot_eq_w(3)-(k3x*w1+k3y*w2);

% matlabFunction(Fb_trans, F3_trans_handle, k3x, k3y, func_c3, 'File', './Functions/AutoFunc_sym_ThrustVerifications', 'Vars', {yaw, alpha, w1, w2});
matlabFunction(Fb_test_w, 'File', './Auto_Functions/AutoFunc_sym_Fb_trans_h', 'Vars', {alpha, Fz, w1, w2});
matlabFunction(F_mot_eq_w(3), 'File', './Auto_Functions/AutoFunc_sym_F3_trans_h', 'Vars', {alpha, Fz, w1, w2});
matlabFunction(k3x, 'File', './Auto_Functions/AutoFunc_sym_k3x_h', 'Vars', {alpha, Fz, w1, w2});
matlabFunction(k3y, 'File', './Auto_Functions/AutoFunc_sym_k3y_h', 'Vars', {alpha, Fz, w1, w2});
matlabFunction(c3, 'File', './Auto_Functions/AutoFunc_sym_c3_h', 'Vars', {alpha, Fz, w1, w2});