
Fd = [0; 0; m_drone_num*g_ref];
alpha_value = 20*pi/180;
Plot_Limits = Lims;

fprop_sat_min = fprop_sat_vec(1);
fprop_sat_max = fprop_sat_vec(2);
N_precision_sphere_cercle = 50;
t = linspace(0, 2*pi, N_precision_sphere_cercle);
if alpha_value < dead_zone_alpha
    alpha_value = 0;
end

[Mixer, ~, Mixer_augmented] = Func_sym_Mixer_wrapper(alpha_value);
pinv_Mixer_augmented = pinv(Mixer_augmented);

Fz_max = (Mixer*fprop_sat_max*ones([8,1])).'*[0; 0; 1; 0; 0; 0];
Fz_min = (Mixer*fprop_sat_min*ones([8,1])).'*[0; 0; 1; 0; 0; 0];
Fz_med = (Fz_max+Fz_min)/2;
Try_vec = linspace(fprop_sat_min, fprop_sat_max, 2);

if Fd(1)~=0 || Fd(2)~=0
    yaw = atan2(Fd(2), Fd(1));
else
    yaw = 0;
end

if yaw == pi/2 || yaw == -pi/2
    Fd_t = Fd(2);
else
    Fd_t = Fd(1)/cos(yaw);
end

if alpha_value < dead_zone_alpha || alpha_value > pi/2 - dead_zone_alpha
    % alpha_value = 0;
    cascade = true;
else
    % alpha_value = alpha_value;
    cascade = false;
end


%% Calcul enveloppe

if ~cascade

    

    P2C_enveloppe = Func_combinaisons(Try_vec, 8);
    n_enveloppe = size(P2C_enveloppe, 2);
    Enveloppe_drone = zeros(6, n_enveloppe);

    for i_e=1:n_enveloppe
        Enveloppe_drone(:, i_e) = Mixer*P2C_enveloppe(:, i_e);
    end

    F_x_enveloppe = Enveloppe_drone(1, :);
    F_y_enveloppe = Enveloppe_drone(2, :);
    F_z_enveloppe = Enveloppe_drone(3, :);

    [tri_enveloppe, ~] = convhull(F_x_enveloppe, F_y_enveloppe, F_z_enveloppe, "Simplify", true);

end


%% Visu sphere Fd

r_sphere = norm(Fd);
[x_sphere, y_sphere, z_sphere] = sphere(N_precision_sphere_cercle);
x_sphere = x_sphere * r_sphere;
y_sphere = y_sphere * r_sphere;
z_sphere = z_sphere * r_sphere;
x_circle_force = r_sphere*cos(t);
y_circle_force = r_sphere*sin(t);

%% Calcul des efforts accessibles
if ~cascade && true

    yaw_eps = pi/20;
    compute_f = (yaw < -yaw_eps && yaw > -pi + yaw_eps) || (yaw > yaw_eps && yaw < pi - yaw_eps);
    compute_c = alpha_value <= pi/4;
    Comb = Func_combinaisons(Try_vec, 2);
    n_comb_v = size(Comb, 2);
    SZ_v = 100;
    P2C_vertical = [Comb, zeros([2, SZ_v-n_comb_v])];
    eps_v = 1e-12;
    i_v = 1;
    compute_v = true;
    while compute_v == true
        % display(P2C_vertical);
        % display(i_v);
        w10_v = P2C_vertical(1,i_v);
        w20_v = P2C_vertical(2,i_v);
        wi0_v = double(AutoFunc_sym_Fmot_trans_wf_c(yaw, w10_v, w20_v, Mixer_augmented, pinv_Mixer_augmented)*compute_f*compute_c + ...
            AutoFunc_sym_Fmot_trans_wl_c(yaw, w10_v, w20_v, Mixer_augmented, pinv_Mixer_augmented)*~compute_f*compute_c + ...
            AutoFunc_sym_Fmot_trans_wf_o(yaw, w10_v, w20_v, Mixer_augmented, pinv_Mixer_augmented)*compute_f*~compute_c + ...
            AutoFunc_sym_Fmot_trans_wl_o(yaw, w10_v, w20_v, Mixer_augmented, pinv_Mixer_augmented)*~compute_f*~compute_c);
        % display(wi0_v);
        
        flag_out_v = false;
        for i_mot = 1:6
            if (i_mot < 1 && compute_c) || (i_mot < 5 && ~compute_c)
                if wi0_v(i_mot) > fprop_sat_max + eps_v
                    prop_sat_v = fprop_sat_max;
                    flag_out_v = true;
                elseif wi0_v(i_mot) < fprop_sat_min - eps_v
                    prop_sat_v = fprop_sat_min;
                    flag_out_v = true;
                else
                    prop_sat_v = 0;
                end
            else
                if wi0_v(i_mot+2) > fprop_sat_max + eps_v
                    prop_sat_v = fprop_sat_max;
                    flag_out_v = true;
                elseif wi0_v(i_mot+2) < fprop_sat_min - eps_v
                    prop_sat_v = fprop_sat_min;
                    flag_out_v = true;
                else
                    prop_sat_v = 0;
                end
            end
            if flag_out_v
                P2C_vertical(:, i_v:end) = [P2C_vertical(:,i_v+1:end), zeros([2,1])];
                % display(P2C_vertical);
                n_comb_v = n_comb_v - 1;
                k1i_v_temp = double(AutoFunc_sym_k1_wf_c(yaw, w10_v, w20_v, Mixer_augmented, pinv_Mixer_augmented)*compute_f*compute_c + ...
                    AutoFunc_sym_k1_wl_c(yaw, w10_v, w20_v, Mixer_augmented, pinv_Mixer_augmented)*~compute_f*compute_c + ...
                    AutoFunc_sym_k1_wf_o(yaw, w10_v, w20_v, Mixer_augmented, pinv_Mixer_augmented)*compute_f*~compute_c + ...
                    AutoFunc_sym_k1_wl_o(yaw, w10_v, w20_v, Mixer_augmented, pinv_Mixer_augmented)*~compute_f*~compute_c);
                k1i_v = double(k1i_v_temp(i_mot));
                k2i_v_temp = double(AutoFunc_sym_k2_wf_c(yaw, w10_v, w20_v, Mixer_augmented, pinv_Mixer_augmented)*compute_f*compute_c + ...
                    AutoFunc_sym_k2_wl_c(yaw, w10_v, w20_v, Mixer_augmented, pinv_Mixer_augmented)*~compute_f*compute_c + ...
                    AutoFunc_sym_k2_wf_o(yaw, w10_v, w20_v, Mixer_augmented, pinv_Mixer_augmented)*compute_f*~compute_c + ...
                    AutoFunc_sym_k2_wl_o(yaw, w10_v, w20_v, Mixer_augmented, pinv_Mixer_augmented)*~compute_f*~compute_c);
                k2i_v = double(k2i_v_temp(i_mot));
                ci_v_temp = double(AutoFunc_sym_c_wf_c(yaw, w10_v, w20_v, Mixer_augmented, pinv_Mixer_augmented)*compute_f*compute_c + ...
                    AutoFunc_sym_c_wl_c(yaw, w10_v, w20_v, Mixer_augmented, pinv_Mixer_augmented)*~compute_f*compute_c + ...
                    AutoFunc_sym_c_wf_o(yaw, w10_v, w20_v, Mixer_augmented, pinv_Mixer_augmented)*compute_f*~compute_c + ...
                    AutoFunc_sym_c_wl_o(yaw, w10_v, w20_v, Mixer_augmented, pinv_Mixer_augmented)*~compute_f*~compute_c);
                ci_v = double(ci_v_temp(i_mot));
                w11_v = w10_v;
                w21_v = (prop_sat_v-k1i_v*w11_v-ci_v)/k2i_v;
                w22_v = w20_v;
                w12_v = (prop_sat_v-k2i_v*w22_v-ci_v)/k1i_v;
                if w21_v >= fprop_sat_min - eps_v && w21_v <= fprop_sat_max + eps_v
                    flag_doublon = false;
                    for i_doublon = 1:n_comb_v
                        if P2C_vertical(1, i_doublon) >= w11_v - eps_v && P2C_vertical(1,i_doublon) <= w11_v + eps_v
                            if P2C_vertical(2, i_doublon) >= w21_v - eps_v && P2C_vertical(2,i_doublon) <= w21_v + eps_v
                                flag_doublon = true;
                            end
                        end
                    end
                    if ~flag_doublon
                        append_v = [w11_v;w21_v];
                        n_comb_v = n_comb_v + 1;
                        P2C_vertical(:, n_comb_v) = append_v;
                        % display(P2C_vertical);
                    end
                end
                if w12_v >= fprop_sat_min - eps_v && w12_v <= fprop_sat_max + eps_v
                    flag_doublon = false;
                    for i_doublon = 1:n_comb_v
                        if P2C_vertical(1, i_doublon) >= w12_v - eps_v && P2C_vertical(1,i_doublon) <= w12_v + eps_v
                            if P2C_vertical(2, i_doublon) >= w22_v - eps_v && P2C_vertical(2,i_doublon) <= w22_v + eps_v
                                flag_doublon = true;
                            end
                        end
                    end
                    if ~flag_doublon
                        append_v = [w12_v; w22_v];
                        n_comb_v = n_comb_v + 1;
                        P2C_vertical(:, n_comb_v) = append_v;
                        % display(P2C_vertical);
                    end
                end
                break;
            end
        end
        if ~flag_out_v
            if i_v < n_comb_v
                i_v = i_v+1;
            else
                compute_v = false;
            end
        end
    end
    P2C_vertical = P2C_vertical(:, 1:n_comb_v);
    n_v = size(P2C_vertical, 2);

    F_trans_computed = zeros(6,n_v);

    for i_v = 1:n_v
        F_trans_computed(:, i_v) = AutoFunc_sym_Fb_trans_wf_c(yaw, P2C_vertical(1,i_v), P2C_vertical(2,i_v), Mixer_augmented, pinv_Mixer_augmented)*compute_f*compute_c + ...
            AutoFunc_sym_Fb_trans_wl_c(yaw, P2C_vertical(1,i_v), P2C_vertical(2,i_v), Mixer_augmented, pinv_Mixer_augmented)*~compute_f*compute_c + ...
            AutoFunc_sym_Fb_trans_wf_o(yaw, P2C_vertical(1,i_v), P2C_vertical(2,i_v), Mixer_augmented, pinv_Mixer_augmented)*compute_f*~compute_c + ...
            AutoFunc_sym_Fb_trans_wl_o(yaw, P2C_vertical(1,i_v), P2C_vertical(2,i_v), Mixer_augmented, pinv_Mixer_augmented)*~compute_f*~compute_c;
    end

    if yaw == pi/2 || yaw == -pi/2
        Ft_v_num = F_trans_computed(2,:);
    else
        Ft_v_num = F_trans_computed(1,:)/cos(yaw);
    end

    Fz_v_num = F_trans_computed(3,:);

    [tri_verticale, ~] = convhull(Ft_v_num, Fz_v_num, "Simplify", true);

    [~, k_bot] = min(Fz_v_num);
    [~, k_top] = max(Fz_v_num);
    [~, k_mid] = max(Ft_v_num);

    %% Calcul si Fd est accessible


    a_side_top = (Fz_v_num(k_mid)-Fz_v_num(k_top))/(Ft_v_num(k_mid)-Ft_v_num(k_top));
    a_side_bot = (Fz_v_num(k_mid)-Fz_v_num(k_bot))/(Ft_v_num(k_mid)-Ft_v_num(k_bot));
    b_side_top = Fz_v_num(k_top);
    b_side_bot = Fz_v_num(k_bot);

    f_side_top = @(t) a_side_top*t + b_side_top;
    f_side_bot = @(t) a_side_bot*t + b_side_bot;
    f_sphere_t = @(t) real((r_sphere.^2-t.^2).^(1/2));

    flag_Fd =  f_side_top(Fd_t)^2+Fd_t^2 >= r_sphere^2 && f_side_bot(Fd_t)^2 + Fd_t^2 <= r_sphere^2;


    %% Calcul des intersections entre le cercle force et la coupe du cone

    if flag_Fd
        Fc_t = Fd_t;
        Fc_z = Fd(3);
        Fc = Fd;
        pitch_c = 0;
        ratio = 1;
    else
        if r_sphere <= Fz_max && r_sphere >= Fz_min
            poly_top = [a_side_top.^2+1, 2*a_side_top*b_side_top, b_side_top.^2-r_sphere.^2];
            poly_bot = [a_side_bot.^2+1, 2*a_side_bot*b_side_bot, b_side_bot.^2-r_sphere.^2];

            roots_top = roots(poly_top);
            roots_bot = roots(poly_bot);
            roots_tot = abs([roots_top; roots_bot]);
            x0 = min(roots_tot);
        else
            x0 = 0;
        end

        if x0 ~= 0
            y0 = f_sphere_t(x0);
        else
            y0 = max([min([Fz_max, r_sphere]), Fz_min]);
        end

        Fc_t = x0;
        Fc_z = y0;
        Fc = [x0*cos(yaw); x0*sin(yaw); y0];
        pitch_c = acos((Fd.'* Fc)/(norm(Fd)*norm(Fc)));
        ratio = norm(Fc)/norm(Fd);
    end
elseif ~cascade && false
    Fc = Fd;
    pitch_c = 0;
    ratio = 1;
else
    x0 = 0;
    y0 = max([min([Fz_max, r_sphere]), Fz_min]);
    Fc_t = x0;
    Fc_z = y0;
    Fc = [x0*cos(yaw); x0*sin(yaw); y0];
    pitch_c = acos((Fd.'* Fc)/(norm(Fd)*norm(Fc)));
    ratio = norm(Fc)/norm(Fd);
end

%% PLOTS

Fx_limits = Plot_Limits(1,:);
Fy_limits = Plot_Limits(2,:);
Fz_limits = Plot_Limits(3,:);
Ft_limits = Plot_Limits(4,:);

frame = figure();
% frame.Visible = 'off';
tiledlayout(1, 2);

nexttile();
hold on
if ~cascade
    trisurf(tri_enveloppe, F_x_enveloppe, F_y_enveloppe, F_z_enveloppe, 'FaceColor', 'c', 'LineStyle', "-", 'EdgeColor', 'k', 'FaceAlpha', 0.2);
end
plot3([0, Fd(1)], [0, Fd(2)], [0, Fd(3)], '-rx');
surf(x_sphere, y_sphere, z_sphere, 'FaceColor', 'r', 'LineStyle', 'none', 'FaceAlpha', 0.3)
plot3([0, Fc(1)], [0, Fc(2)], [0, Fc(3)], '--go')
hold off
axis equal
grid();
xlabel('$F_x$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$F_y$', 'Interpreter', 'latex', 'FontSize', 20);
zlabel('$F_z$', 'Interpreter', 'latex', 'FontSize', 20, 'Rotation', 0);
title("Volume des forces atteignables pour $\alpha=$"+num2str(alpha_value*180/pi)+"$^\circ$", 'Interpreter', 'latex')
if ~cascade
    legend("Enveloppe forces atteignables", ...
        "$F_d$", "Sphere des rotations possibles de $F_d$", "Rotation optimale de $F_d$ de "+num2str(pitch_c*180/pi)+"°", 'Interpreter', 'latex', 'Location', 'northeast')
else
    legend("$F_d$", "Sphere des rotations possibles de $F_d$", "Rotation optimale de $F_d$ de "+num2str(pitch_c*180/pi)+"°", 'Interpreter', 'latex', 'Location', 'northeast')
end
view(45, 30)
lightangle(-45,30)
lighting flat
xlim(Fx_limits);
ylim(Fy_limits);
zlim(Fz_limits);

nexttile();
hold on
if ~cascade
    plot(Ft_v_num(tri_verticale), Fz_v_num(tri_verticale), "c");
end
plot([0, Fd_t], [0, Fd(3)], '-rx')
plot(x_circle_force, y_circle_force, '--r')
plot([0, Fc_t], [0, Fc_z], '--go')
hold off
axis equal
grid()
xlabel('$F_t$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$F_z$', 'Interpreter', 'latex', 'FontSize', 20);
title("Volume des forces atteignables pour $\alpha=$"+num2str(alpha_value*180/pi)+"$^\circ$, dans le plan vertical de l'effort horizontal", 'Interpreter', 'latex')
if ~cascade
    legend("Forces atteignables pour $\psi="+num2str(yaw*180/pi)+"$°", "$F_d$", "Cercle des rotations possibles de $F_d$", ...
        "Rotation optimale de $F_d$ de "+num2str(pitch_c*180/pi)+"°", 'Interpreter', 'latex', 'Location', 'northeast')
else
    legend("$F_d$", "Cercle des rotations possibles de $F_d$", ...
        "Rotation optimale de $F_d$ de "+num2str(pitch_c*180/pi)+"°", 'Interpreter', 'latex', 'Location', 'northeast')
end
xlim(Ft_limits);
ylim(Fz_limits);