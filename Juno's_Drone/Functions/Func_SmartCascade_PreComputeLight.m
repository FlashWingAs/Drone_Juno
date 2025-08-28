function [kFmot_lambda] = Func_SmartCascade_PreComputeLight(alpha, yaw)

%% Calculs symboliques

% text = 'Pré-calculs symboliques du controlleur';
% Func_Progress_Log(false, text, LOG_SIZE);

% Espace de forces moteurs pour un torseur donné
syms Ft Fz
n_actionneur = 8;

Fx = cos(yaw)*Ft;
Fy = sin(yaw)*Ft;
Fb = [Fx; Fy; Fz; 0; 0; 0];
J_mot = AutoFunc_sym_Mixer(alpha);
F_mot = pinv(J_mot)*Fb;

% Lien des efforts tangeantiels par les vitesses moteurs

syms lambda_1 lambda_2

eps_equal = 1e-1;
sys_found = false;
i = 1;
k_Ft = -1;
k_Fz = -1;
while ~sys_found
    if k_Fz < 0
        var_Fz = double(diff(F_mot(i), Fz));
        if ~Func_equal_eps(var_Fz, 0, eps_equal)
            k_Fz = i;
            i = 1;
        end
    end
    if k_Fz > 0 && i ~= k_Fz
        var_Ft = double(diff(F_mot(i), Ft));
        if ~Func_equal_eps(var_Ft, 0, eps_equal)
            k_Ft = i;
            sys_found = true;
        end
    end
    if i >= n_actionneur
        i = 1;
    else
        i = i+1;
    end
end

F_mot_Ft = F_mot(k_Ft);
F_mot_Fz = F_mot(k_Fz);
Fz_lambda = solve(F_mot_Fz == lambda_1, Fz);
F_mot_Ft = subs(F_mot_Ft, Fz, Fz_lambda);
Ft_lambda = solve(F_mot_Ft == lambda_2, Ft);
Fz_lambda = subs(Fz_lambda, Ft, Ft_lambda);

F_mot_lambda = simplify(subs(F_mot, [Fz, Ft], [Fz_lambda, Ft_lambda]));

kFmot_lambda = zeros([2,n_actionneur]);
for i = 1:n_actionneur
    kFmot_lambda(1, i) = double(diff(F_mot_lambda(i), lambda_1));
    kFmot_lambda(2, i) = double(diff(F_mot_lambda(i), lambda_2));
end


end