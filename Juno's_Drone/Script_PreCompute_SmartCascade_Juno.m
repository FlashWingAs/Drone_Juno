%% Calculs symboliques

text = 'Pré-calculs symboliques du controlleur';
Func_Progress_Log(false, text, LOG_SIZE);

syms Ft Fz yaw_sym
Mixer_augmented = sym('Mixer_augmented_%d%d', [8, 8]);
pinv_Mixer_augmented = sym('pinv_Mixer_augmented_%d%d', [8, 8]);
k1f_c = sym('k1%df', [6,1]);
k2f_c = sym('k2%df', [6,1]);
cf_c = sym('c%df', [6,1]);
k1l_c = sym('k1%dl', [6,1]);
k2l_c = sym('k2%dl', [6,1]);
cl_c = sym('c%dl', [6,1]);
k1f_o = sym('k1%df', [6,1]);
k2f_o = sym('k2%df', [6,1]);
cf_o = sym('c%df', [6,1]);
k1l_o = sym('k1%dl', [6,1]);
k2l_o = sym('k2%dl', [6,1]);
cl_o = sym('c%dl', [6,1]);

Fx = cos(yaw_sym)*Ft;
Fy = sin(yaw_sym)*Ft;
Fb_test = [Fx; Fy; Fz; 0; 0; 0; 0; 0];
F_mot = pinv_Mixer_augmented*Fb_test;
F_mot_eq = formula(F_mot);

Fmot1 = F_mot_eq(1);
Fmot2 = F_mot_eq(2);

Fmot5 = F_mot_eq(5);
Fmot6 = F_mot_eq(6);

% Lien des efforts tangeantiels par les vitesses moteurs

syms w1 w2

yaw_eps = pi/20;

% Première approche

Fz_wf_c = solve(Fmot2==w1, Fz);
Fmot1_wf_c = subs(Fmot1, Fz, Fz_wf_c);
Ft_wf_c = solve(Fmot1_wf_c==w2, Ft);
F_mot_eq_wf_c = simplify(subs(subs(F_mot_eq, Fz, Fz_wf_c), Ft, Ft_wf_c));

Fz_wl_c = solve(Fmot1==w1, Fz);
Fmot2_wl_c = subs(Fmot2, Fz, Fz_wl_c);
Ft_wl_c = solve(Fmot2_wl_c==w2, Ft);
F_mot_eq_wl_c = simplify(subs(subs(F_mot_eq, Fz, Fz_wl_c), Ft, Ft_wl_c));

Fb_test_wf_c = Mixer_augmented*F_mot_eq_wf_c;
Fb_test_wl_c = Mixer_augmented*F_mot_eq_wl_c;

func_k1f_c = cell(6,1);
func_k2f_c = cell(6,1);
func_cf_c = cell(6,1);
func_k1l_c = cell(6,1);
func_k2l_c = cell(6,1);
func_cl_c = cell(6,1);

for i=1:6
    k1f_c(i) = diff(F_mot_eq_wf_c(i+2), w1);
    k2f_c(i) = diff(F_mot_eq_wf_c(i+2), w2);
    cf_c(i) = F_mot_eq_wf_c(i+2)-(k1f_c(i)*w1+k2f_c(i)*w2);
    k1l_c(i) = diff(F_mot_eq_wl_c(i+2), w1);
    k2l_c(i) = diff(F_mot_eq_wl_c(i+2), w2);
    cl_c(i) = F_mot_eq_wl_c(i+2)-(k1l_c(i)*w1+k2l_c(i)*w2);
end

% Seconde approche

Fz_wf_o = solve(Fmot6==w1, Fz);
Fmot1_wf_o = subs(Fmot5, Fz, Fz_wf_o);
Ft_wf_o = solve(Fmot1_wf_o==w2, Ft);
F_mot_eq_wf_o = simplify(subs(subs(F_mot_eq, Fz, Fz_wf_o), Ft, Ft_wf_o));

Fz_wl_o = solve(Fmot5==w1, Fz);
Fmot2_wl_o = subs(Fmot6, Fz, Fz_wl_o);
Ft_wl_o = solve(Fmot2_wl_o==w2, Ft);
F_mot_eq_wl_o = simplify(subs(subs(F_mot_eq, Fz, Fz_wl_o), Ft, Ft_wl_o));

Fb_test_wf_o = Mixer_augmented*F_mot_eq_wf_o;
Fb_test_wl_o = Mixer_augmented*F_mot_eq_wl_o;

func_k1f_o = cell(6,1);
func_k2f_o = cell(6,1);
func_cf_o = cell(6,1);
func_k1l_o = cell(6,1);
func_k2l_o = cell(6,1);
func_cl_o = cell(6,1);

for i=1:6
    if i<5
        k1f_o(i) = diff(F_mot_eq_wf_o(i), w1);
        k2f_o(i) = diff(F_mot_eq_wf_o(i), w2);
        cf_o(i) = F_mot_eq_wf_o(i)-(k1f_o(i)*w1+k2f_o(i)*w2);
        k1l_o(i) = diff(F_mot_eq_wl_o(i), w1);
        k2l_o(i) = diff(F_mot_eq_wl_o(i), w2);
        cl_o(i) = F_mot_eq_wl_o(i)-(k1l_o(i)*w1+k2l_o(i)*w2);
    else
        k1f_o(i) = diff(F_mot_eq_wf_o(i+2), w1);
        k2f_o(i) = diff(F_mot_eq_wf_o(i+2), w2);
        cf_o(i) = F_mot_eq_wf_o(i+2)-(k1f_o(i)*w1+k2f_o(i)*w2);
        k1l_o(i) = diff(F_mot_eq_wl_o(i+2), w1);
        k2l_o(i) = diff(F_mot_eq_wl_o(i+2), w2);
        cl_o(i) = F_mot_eq_wl_o(i+2)-(k1l_o(i)*w1+k2l_o(i)*w2);
    end
end

Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

%% Génération des fonctions

text = 'Génération des fonctions internes du controlleur';
Func_Progress_Log(false, text, LOG_SIZE);

matlabFunction(Fb_test_wf_c(1:6), 'File', './Auto_Functions/AutoFunc_sym_Fb_trans_wf_c', 'Vars', {yaw_sym, w1, w2, Mixer_augmented, pinv_Mixer_augmented});
matlabFunction(F_mot_eq_wf_c, 'File', './Auto_Functions/AutoFunc_sym_Fmot_trans_wf_c', 'Vars', {yaw_sym, w1, w2, Mixer_augmented, pinv_Mixer_augmented});
matlabFunction(k1f_c, 'File', './Auto_Functions/AutoFunc_sym_k1_wf_c', 'Vars', {yaw_sym, w1, w2, Mixer_augmented, pinv_Mixer_augmented});
matlabFunction(k2f_c, 'File', './Auto_Functions/AutoFunc_sym_k2_wf_c', 'Vars', {yaw_sym, w1, w2, Mixer_augmented, pinv_Mixer_augmented});
matlabFunction(cf_c, 'File', './Auto_Functions/AutoFunc_sym_c_wf_c', 'Vars', {yaw_sym, w1, w2, Mixer_augmented, pinv_Mixer_augmented});

matlabFunction(Fb_test_wl_c(1:6), 'File', './Auto_Functions/AutoFunc_sym_Fb_trans_wl_c', 'Vars', {yaw_sym, w1, w2, Mixer_augmented, pinv_Mixer_augmented});
matlabFunction(F_mot_eq_wl_c, 'File', './Auto_Functions/AutoFunc_sym_Fmot_trans_wl_c', 'Vars', {yaw_sym, w1, w2, Mixer_augmented, pinv_Mixer_augmented});
matlabFunction(k1l_c, 'File', './Auto_Functions/AutoFunc_sym_k1_wl_c', 'Vars', {yaw_sym, w1, w2, Mixer_augmented, pinv_Mixer_augmented});
matlabFunction(k2l_c, 'File', './Auto_Functions/AutoFunc_sym_k2_wl_c', 'Vars', {yaw_sym, w1, w2, Mixer_augmented, pinv_Mixer_augmented});
matlabFunction(cl_c, 'File', './Auto_Functions/AutoFunc_sym_c_wl_c', 'Vars', {yaw_sym, w1, w2, Mixer_augmented, pinv_Mixer_augmented});


matlabFunction(Fb_test_wf_o(1:6), 'File', './Auto_Functions/AutoFunc_sym_Fb_trans_wf_o', 'Vars', {yaw_sym, w1, w2, Mixer_augmented, pinv_Mixer_augmented});
matlabFunction(F_mot_eq_wf_o, 'File', './Auto_Functions/AutoFunc_sym_Fmot_trans_wf_o', 'Vars', {yaw_sym, w1, w2, Mixer_augmented, pinv_Mixer_augmented});
matlabFunction(k1f_o, 'File', './Auto_Functions/AutoFunc_sym_k1_wf_o', 'Vars', {yaw_sym, w1, w2, Mixer_augmented, pinv_Mixer_augmented});
matlabFunction(k2f_o, 'File', './Auto_Functions/AutoFunc_sym_k2_wf_o', 'Vars', {yaw_sym, w1, w2, Mixer_augmented, pinv_Mixer_augmented});
matlabFunction(cf_o, 'File', './Auto_Functions/AutoFunc_sym_c_wf_o', 'Vars', {yaw_sym, w1, w2, Mixer_augmented, pinv_Mixer_augmented});

matlabFunction(Fb_test_wl_o(1:6), 'File', './Auto_Functions/AutoFunc_sym_Fb_trans_wl_o', 'Vars', {yaw_sym, w1, w2, Mixer_augmented, pinv_Mixer_augmented});
matlabFunction(F_mot_eq_wl_o, 'File', './Auto_Functions/AutoFunc_sym_Fmot_trans_wl_o', 'Vars', {yaw_sym, w1, w2, Mixer_augmented, pinv_Mixer_augmented});
matlabFunction(k1l_o, 'File', './Auto_Functions/AutoFunc_sym_k1_wl_o', 'Vars', {yaw_sym, w1, w2, Mixer_augmented, pinv_Mixer_augmented});
matlabFunction(k2l_o, 'File', './Auto_Functions/AutoFunc_sym_k2_wl_o', 'Vars', {yaw_sym, w1, w2, Mixer_augmented, pinv_Mixer_augmented});
matlabFunction(cl_o, 'File', './Auto_Functions/AutoFunc_sym_c_wl_o', 'Vars', {yaw_sym, w1, w2, Mixer_augmented, pinv_Mixer_augmented});

Func_Progress_Log(true, text, LOG_SIZE);
current_step = current_step + 1;
waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")