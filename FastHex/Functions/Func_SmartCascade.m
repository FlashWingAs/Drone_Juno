function [Fc, Rot_Mat, ratio] = Func_SmartCascade_V3(Fd, Mixer, alpha, sat_vec, dead_zone_alpha)

sat_max = sat_vec(2);
sat_min = sat_vec(1);
yaw = atan2(Fd(2), Fd(1));
if yaw == pi/2
    Fd_t = Fd(2);
elseif yaw == -pi/2
    Fd_t = -Fd(2);
else
    Fd_t = Fd(1)/cos(yaw);
end
Fz_max = (Mixer*sat_max*ones([6,1])).'*[0; 0; 1; 0; 0; 0];
Fz_min = (Mixer*sat_min*ones([6,1])).'*[0; 0; 1; 0; 0; 0];
Fz_med = (Fz_max+Fz_min)/2;
r_sphere = norm(Fd);


if alpha < dead_zone_alpha
    alpha_value = 0;
    cascade = true;
else
    alpha_value = alpha;
    cascade = false;
end

%% Calcul du cone
if ~cascade
    Try_vec = sat_vec;
    Comb = Func_combinaisons(Try_vec, 2);
    n_comb = size(Comb, 2);
    SZ = 20;
    Points2compute = [Comb, zeros([2, SZ-n_comb])];
    eps = 1e-12;
    i = 1;
    compute = true;
    eps_yaw = pi/20;
    compute_via_w3 =  (yaw > -pi/3 - eps_yaw && yaw < -pi/3 + eps_yaw) || (yaw > 2*pi/3 - eps_yaw && yaw < 2*pi/3 + eps_yaw);
    while compute == true
        x0 = Points2compute(1,i);
        y0 = Points2compute(2,i);
        z0 = AutoFunc_sym_F3_trans_w2(alpha_value, yaw, x0, y0)*~compute_via_w3 + AutoFunc_sym_F3_trans_w3(alpha_value, yaw, x0, y0)*compute_via_w3;
        flag_out = false;
        if z0 > sat_max + eps
            prop_sat = sat_max;
            flag_out = true;
        elseif z0 < sat_min - eps
            prop_sat = sat_min;
            flag_out = true;
        else
            prop_sat = 0;
        end
        if flag_out
            Points2compute(:, i:end) = [Points2compute(:,i+1:end), zeros([2,1])];
            n_comb = n_comb - 1;
            k3x = AutoFunc_sym_k3x_w2(alpha_value, yaw, x0, y0)*~compute_via_w3 + AutoFunc_sym_k3x_w3(alpha_value, yaw, x0, y0)*compute_via_w3;
            k3y = AutoFunc_sym_k3y_w2(alpha_value, yaw, x0, y0)*~compute_via_w3 + AutoFunc_sym_k3y_w3(alpha_value, yaw, x0, y0)*compute_via_w3;
            c3 = AutoFunc_sym_c3_w2(alpha_value, yaw, x0, y0)*~compute_via_w3 + AutoFunc_sym_c3_w3(alpha_value, yaw, x0, y0)*compute_via_w3;
            y1 = (prop_sat - c3 - k3x*x0)/k3y;
            x1 = (prop_sat - c3 - k3y*y0)/k3x;
            z10 = AutoFunc_sym_F3_trans_w2(alpha_value, yaw, x1, y0)*~compute_via_w3 + AutoFunc_sym_F3_trans_w3(alpha_value, yaw, x1, y0)*compute_via_w3;
            z01 = AutoFunc_sym_F3_trans_w2(alpha_value, yaw, x0, y1)*~compute_via_w3 + AutoFunc_sym_F3_trans_w3(alpha_value, yaw, x0, y1)*compute_via_w3;
            if x1 >= sat_min - eps && x1 <= sat_max + eps
                if z10 >= sat_min-eps && z10 <= sat_max + eps
                    append = [x1;y0];
                    n_comb = n_comb + 1;
                    Points2compute(:,n_comb) = append;
                end
            end
            if y1 >= sat_min - eps && y1 <= sat_max + eps
                if z01 >= sat_min - eps && z01 <= sat_max + eps
                    append = [x0; y1];
                    n_comb = n_comb + 1;
                    Points2compute(:,n_comb) = append;
                end
            end
            if i > n_comb
                compute = false;
            end
        else
            if i < n_comb
                i = i+1;
            else
                compute = false;
            end
        end
    end
    PointsComputed = Points2compute(:,1:n_comb);
    n = size(PointsComputed, 2);

    if n > 0 && size(find(PointsComputed), 1) > 0
        F_trans_computed = zeros(6,n);
        F3_test_w = zeros([1,n]);
        for i = 1:n
            F_trans_computed(:, i) = AutoFunc_sym_Fb_trans_w2(alpha_value, yaw, PointsComputed(1,i), PointsComputed(2,i))*~compute_via_w3 + AutoFunc_sym_Fb_trans_w3(alpha_value, yaw, PointsComputed(1,i), PointsComputed(2,i))*compute_via_w3;
            F3_test_w(i) = AutoFunc_sym_F3_trans_w2(alpha_value, yaw, PointsComputed(1,i), PointsComputed(2,i))*~compute_via_w3 + AutoFunc_sym_F3_trans_w3(alpha_value, yaw, PointsComputed(1,i), PointsComputed(2,i))*compute_via_w3;
        end
        if yaw == pi/2 || yaw == -pi/2
            Ft_num = F_trans_computed(2,:);
        else    
            Ft_num = F_trans_computed(1,:)/cos(yaw);
        end
        Fz_num = F_trans_computed(3,:);
    else
        Ft_num = zeros(1,1);
        Fz_num = zeros(1,1);
    end

    [~, k_bot] = min(Fz_num);
    [~, k_top] = max(Fz_num);
    [~, k_mid] = max(Ft_num);

    %% Calcul si Fd est accessible

    A = [Ft_num(k_top); Fz_num(k_top)];
    B = [Ft_num(k_mid); Fz_num(k_mid)];
    C = [Ft_num(k_bot); Fz_num(k_bot)];

    a_side_top = (B(2)-A(2))/(B(1)-A(1));
    a_side_bot = (B(2)-C(2))/(B(1)-C(1));
    b_side_top = A(2);
    b_side_bot = C(2);

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