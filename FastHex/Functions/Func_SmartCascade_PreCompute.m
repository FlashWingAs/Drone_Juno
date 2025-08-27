function [F_mot_lambda, Fb_lambda, kFmot_lambda, kFb_lambda] = Func_SmartCascade_PreCompute(alpha, yaw, n_actionneur)

%% Calculs symboliques

% text = 'Pré-calculs symboliques du controlleur';
% Func_Progress_Log(false, text, LOG_SIZE);

% Espace de forces moteurs pour un torseur donné
syms Ft Fz

Fx = cos(yaw)*Ft;
Fy = sin(yaw)*Ft;
Fb = [Fx; Fy; Fz; 0; 0; 0];
J_mot = AutoFunc_sym_Mixer(alpha);
F_mot = J_mot\Fb;

% Lien des efforts tangeantiels par les vitesses moteurs

syms lambda_1 lambda_2

eps_equal = 1e-5;
sys_found = false;
i = 1;
k_Ft = -1;
k_Fz = -1;
while ~sys_found
    if k_Ft < 0
        var_Ft = double(diff(F_mot(i), Ft));
        if ~Func_equal_eps(var_Ft, 0, eps_equal)
            k_Ft = i;
            i = 1;
        end
    end
    if k_Ft > 0 && i ~= k_Ft
        var_Fz = double(diff(F_mot(i), Fz));
        if ~Func_equal_eps(var_Fz, 0, eps_equal)
            k_Fz = i;
            sys_found = true;
        end
    end
    if i >= n_actionneur
        i = 1;
    else
        i = i+1;
    end
end

F_mot_Fz = F_mot(k_Fz);
F_mot_Ft = F_mot(k_Ft);
Ft_lambda = solve(F_mot_Ft == lambda_1, Ft);
F_mot_Fz = subs(F_mot_Fz, Ft, Ft_lambda);
Fz_lambda = solve(F_mot_Fz == lambda_2, Fz);
Ft_lambda = subs(Ft_lambda, Fz, Fz_lambda);

F_mot_lambda = simplify(subs(F_mot, [Ft, Fz], [Ft_lambda, Fz_lambda]));
Fb_lambda = J_mot*F_mot_lambda;

kFmot_lambda = cell([2,n_actionneur]);
kFb_lambda = cell([2,6]);
for i = 1:n_actionneur
    kFmot_lambda{1, i} = double(diff(F_mot_lambda(i), lambda_1));
    kFmot_lambda{2, i} = double(diff(F_mot_lambda(i), lambda_2));
end
for i = 1:6
    kFb_lambda{1, i} = double(diff(Fb_lambda(i), lambda_1));
    kFb_lambda{2, i} = double(diff(Fb_lambda(i), lambda_2));
end


% Func_Progress_Log(true, text, LOG_SIZE);
% current_step = current_step + 1;
% waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")

%% Génération des fonctions

% text = 'Génération des fonctions internes du controlleur';
% Func_Progress_Log(false, text, LOG_SIZE);

% matlabFunction(kFmot_lambda, 'File', './Auto_Functions/AutoFunc_kFmot_lambda', 'Vars', {lambda_1, lambda_2});
% matlabFunction(kFb_lambda, 'File', './Auto_Functions/AutoFunc_kFb_lambda', 'Vars', {lambda_1, lambda_2});
% matlabFunction(F_mot_lambda, 'File', './Auto_Functions/AutoFunc_F_mot_lambda', 'Vars', {lambda_1, lambda_2});
% matlabFunction(kFb_lambda, 'File', './Auto_Functions/AutoFunc_kFb_lambda', 'Vars', {lambda_1, lambda_2});

% Func_Progress_Log(true, text, LOG_SIZE);
% current_step = current_step + 1;
% waitbar(current_step/number_of_steps, objWaitBar, "calcul en cours - "+num2str(current_step/number_of_steps*100,3)+"%")