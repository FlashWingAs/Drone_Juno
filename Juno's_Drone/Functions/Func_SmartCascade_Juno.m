function [Fc, Rot_Mat, ratio] = Func_SmartCascade_Juno(Fd, Mixer, Mixer_augmented, pinv_Mixer_augmented, alpha, sat_vec, dead_zone_alpha)

fprop_sat_min = sat_vec(1);
fprop_sat_max = sat_vec(2);
yaw = atan2(Fd(2), Fd(1));

if yaw == pi/2
    Fd_t = Fd(2);
elseif yaw == -pi/2
    Fd_t = -Fd(2);
else
    Fd_t = Fd(1)/cos(yaw);
end
Fz_max = (Mixer*fprop_sat_max*ones([8,1])).'*[0; 0; 1; 0; 0; 0];
Fz_min = (Mixer*fprop_sat_min*ones([8,1])).'*[0; 0; 1; 0; 0; 0];
Fz_med = (Fz_max+Fz_min)/2;
r_sphere = norm(Fd);
Try_vec = sat_vec;


if alpha < dead_zone_alpha || alpha > pi/2 - dead_zone_alpha
    % alpha_value = 0;
    cascade = true;
else
    % alpha_value = alpha;
    cascade = false;
end

%% Calcul des efforts accessibles
if ~cascade && true

    yaw_eps = pi/20;
    compute_f = (yaw < -yaw_eps && yaw > -pi + yaw_eps) || (yaw > yaw_eps && yaw < pi - yaw_eps);
    compute_c = alpha <= pi/4;
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
        Ft_num = F_trans_computed(2,:);
    else
        Ft_num = F_trans_computed(1,:)/cos(yaw);
    end

    Fz_num = F_trans_computed(3,:);

    [~, k_bot] = min(Fz_num);
    [~, k_top] = max(Fz_num);
    [~, k_mid] = max(Ft_num);

    %% Calcul si Fd est accessible


    a_side_top = (Fz_num(k_mid)-Fz_num(k_top))/(Ft_num(k_mid)-Ft_num(k_top));
    a_side_bot = (Fz_num(k_mid)-Fz_num(k_bot))/(Ft_num(k_mid)-Ft_num(k_bot));
    b_side_top = Fz_num(k_top);
    b_side_bot = Fz_num(k_bot);

    f_side_top = @(t) a_side_top*t + b_side_top;
    f_side_bot = @(t) a_side_bot*t + b_side_bot;
    f_sphere_t = @(t) real((r_sphere.^2-t.^2).^(1/2));

    flag_Fd =  f_side_top(Fd_t)^2+Fd_t^2 >= r_sphere^2 && f_side_bot(Fd_t)^2 + Fd_t^2 <= r_sphere^2;


    %% Calcul des intersections entre le cercle force et la coupe du cone

    if flag_Fd
        % Fc_t = Fd_t;
        % Fc_z = Fd(3);
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

        % Fc_t = x0;
        % Fc_z = y0;
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
    % Fc_t = x0;
    % Fc_z = y0;
    Fc = [x0*cos(yaw); x0*sin(yaw); y0];
    pitch_c = acos((Fd.'* Fc)/(norm(Fd)*norm(Fc)));
    ratio = norm(Fc)/norm(Fd);
end

Rot_Mat = Func_rot_z(yaw)*Func_rot_y(pitch_c)*Func_rot_z(-yaw);

end