function [frame] = Func_visu_movie(Fd, Mixer_flat, alpha_value, dead_zone_alpha, fprop_sat_vec, Plot_Limits)
fprop_sat_min = fprop_sat_vec(1);
fprop_sat_max = fprop_sat_vec(2);
N_precision_sphere_cercle = 50;
t = linspace(0, 2*pi, N_precision_sphere_cercle);
if alpha_value < dead_zone_alpha
    alpha_value = 0;
end

Mixer = Func_sym_Mixer_wrapper(alpha_value, dead_zone_alpha, Mixer_flat);

Fz_max = (Mixer*fprop_sat_max*ones([6,1])).'*[0; 0; 1; 0; 0; 0];
Fz_min = (Mixer*fprop_sat_min*ones([6,1])).'*[0; 0; 1; 0; 0; 0];
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

if alpha_value < dead_zone_alpha
    alpha_value = 0;
    cascade = true;
else
    cascade = false;
end

%% Calcul enveloppe

if ~cascade

    P2C_enveloppe = Func_combinaisons(Try_vec, 3);
    n_enveloppe = size(P2C_enveloppe, 2);
    Enveloppe_drone = zeros(6, n_enveloppe);

    for i_e=1:n_enveloppe
        Enveloppe_drone(:, i_e) = Mixer*[P2C_enveloppe(:, i_e); P2C_enveloppe(:, i_e)];
    end

    F_x_enveloppe = Enveloppe_drone(1, :);
    F_y_enveloppe = Enveloppe_drone(2, :);
    F_z_enveloppe = Enveloppe_drone(3, :);

    [tri_enveloppe, ~] = convhull(F_x_enveloppe, F_y_enveloppe, F_z_enveloppe, "Simplify", true);

    Comb = Func_combinaisons(Try_vec, 2);
    n_comb_h = size(Comb, 2);
    SZ_h = 20;
    P2C_horizontal = [Comb, zeros([2, SZ_h-n_comb_h])];
    eps_h = 1e-12;
    i_h = 1;
    compute_h = true;
    while compute_h == true
        x0_h = P2C_horizontal(1,i_h);
        y0_h = P2C_horizontal(2,i_h);
        z0_h = AutoFunc_sym_F3_trans_h(alpha_value, Fz_med, x0_h, y0_h);
        flag_out_h = false;
        if z0_h > fprop_sat_max+ eps_h
            prop_sat_h = fprop_sat_max;
            flag_out_h = true;
        elseif z0_h < fprop_sat_min - eps_h
            prop_sat_h = fprop_sat_min;
            flag_out_h = true;
        end
        if flag_out_h
            P2C_horizontal(:, i_h:end) = [P2C_horizontal(:,i_h+1:end), zeros([2,1])];
            n_comb_h = n_comb_h - 1;
            k3x_h = AutoFunc_sym_k3x_h(alpha_value, Fz_med, x0_h, y0_h);
            k3y_h = AutoFunc_sym_k3y_h(alpha_value, Fz_med, x0_h, y0_h);
            c3_h = AutoFunc_sym_c3_h(alpha_value, Fz_med, x0_h, y0_h);
            y1_h = (prop_sat_h - c3_h - k3x_h*x0_h)/k3y_h;
            x1_h = (prop_sat_h - c3_h - k3y_h*y0_h)/k3x_h;
            z10_h = AutoFunc_sym_F3_trans_h(alpha_value, Fz_med, x1_h, y0_h);
            z01_h = AutoFunc_sym_F3_trans_h(alpha_value, Fz_med, x0_h, y1_h);
            if x1_h >= fprop_sat_min - eps_h && x1_h <= fprop_sat_max + eps_h
                if z10_h >= fprop_sat_min - eps_h && z10_h <= fprop_sat_max + eps_h
                    append_h = [x1_h;y0_h];
                    n_comb_h = n_comb_h + 1;
                    P2C_horizontal(:, n_comb_h) = append_h;
                end
            end
            if y1_h >= fprop_sat_min - eps_h && y1_h <= fprop_sat_max + eps_h
                if z01_h >= fprop_sat_min - eps_h && z01_h <= fprop_sat_max + eps_h
                    append_h = [x0_h; y1_h];
                    n_comb_h = n_comb_h + 1;
                    P2C_horizontal(:, n_comb_h) = append_h;
                end
            end
            if i_h > n_comb_h
                compute_h = false;
            end
        else
            if i_h < n_comb_h
                i_h = i_h+1;
            else
                compute_h = false;
            end
        end
    end
    P2C_horizontal = P2C_horizontal(:, 1:n_comb_h);
    n_h = size(P2C_horizontal, 2);

    if n_h > 0 && size(find(P2C_horizontal), 1) > 0

        F_trans_computed_h = zeros(6,n_h);
        % F3_test_h_w = zeros([1,n_h]);

        for i_h = 1:n_h
            F_trans_computed_h(:, i_h) = AutoFunc_sym_Fb_trans_h(alpha_value, Fz_med, P2C_horizontal(1,i_h), P2C_horizontal(2,i_h));
            % F3_test_h_w(i_h) = AutoFunc_sym_F3_trans_h(alpha_value, Fz_med, P2C_horizontal(1,i_h), P2C_horizontal(2,i_h));
        end

        Fy_h_num = F_trans_computed_h(2,:);
        a_hexa = max(Fy_h_num);
        r_cone = sqrt(3)/2*a_hexa;
        [x_cone, y_cone, z_cone] = cylinder([0, r_cone, 0], 50);
        z_cone = z_cone*(max(F_z_enveloppe)-min(F_z_enveloppe));

    end

end


%% Visu sphere Fd

r_sphere = norm(Fd);
[x_sphere, y_sphere, z_sphere] = sphere(N_precision_sphere_cercle);
x_sphere = x_sphere * r_sphere;
y_sphere = y_sphere * r_sphere;
z_sphere = z_sphere * r_sphere;
x_circle_force = r_sphere*cos(t);
y_circle_force = r_sphere*sin(t);

%% Calcul slice verticale plan Fd


if ~cascade

    n_comb_v = size(Comb, 2);
    SZ_v = 20;
    P2C_vertical = [Comb, zeros([2, SZ_v-n_comb_v])];
    eps_v = 1e-12;
    i_v = 1;
    compute_v = true;
    eps_yaw = pi/20;
    compute_via_w3 =  (yaw > -pi/3 - eps_yaw && yaw < -pi/3 + eps_yaw) || (yaw > 2*pi/3 - eps_yaw && yaw < 2*pi/3 + eps_yaw);
    while compute_v == true
        x0_v = P2C_vertical(1,i_v);
        y0_v = P2C_vertical(2,i_v);
        z0_v = AutoFunc_sym_F3_trans_w2(alpha_value, yaw, x0_v, y0_v)*~compute_via_w3 + AutoFunc_sym_F3_trans_w3(alpha_value, yaw, x0_v, y0_v)*compute_via_w3;
        flag_out_v = false;
        if z0_v > fprop_sat_max+ eps_v
            prop_sat_v = fprop_sat_max;
            flag_out_v = true;
        elseif z0_v < fprop_sat_min - eps_h
            prop_sat_v = fprop_sat_min;
            flag_out_v = true;
        end
        if flag_out_v
            P2C_vertical(:, i_v:end) = [P2C_vertical(:,i_v+1:end), zeros([2,1])];
            n_comb_v = n_comb_v - 1;
            k3x_v = AutoFunc_sym_k3x_w2(alpha_value, yaw, x0_v, y0_v)*~compute_via_w3 + AutoFunc_sym_k3x_w3(alpha_value, yaw, x0_v, y0_v)*compute_via_w3;
            k3y_v = AutoFunc_sym_k3y_w2(alpha_value, yaw, x0_v, y0_v)*~compute_via_w3 + AutoFunc_sym_k3y_w3(alpha_value, yaw, x0_v, y0_v)*compute_via_w3;
            c3_v = AutoFunc_sym_c3_w2(alpha_value, yaw, x0_v, y0_v)*~compute_via_w3 + AutoFunc_sym_c3_w3(alpha_value, yaw, x0_v, y0_v)*compute_via_w3;
            y1_v = (prop_sat_v - c3_v - k3x_v*x0_v)/k3y_v;
            x1_v = (prop_sat_v - c3_v - k3y_v*y0_v)/k3x_v;
            z10_v = AutoFunc_sym_F3_trans_w2(alpha_value, yaw, x1_v, y0_v)*~compute_via_w3 + AutoFunc_sym_F3_trans_w3(alpha_value, yaw, x1_v, y0_v)*compute_via_w3;
            z01_v = AutoFunc_sym_F3_trans_w2(alpha_value, yaw, x0_v, y1_v)*~compute_via_w3 + AutoFunc_sym_F3_trans_w3(alpha_value, yaw, x0_v, y1_v)*compute_via_w3;
            if x1_v >= fprop_sat_min - eps_v && x1_v <= fprop_sat_max + eps_v
                if z10_v >= fprop_sat_min - eps_v && z10_v <= fprop_sat_max + eps_v
                    append_v = [x1_v;y0_v];
                    n_comb_v = n_comb_v + 1;
                    P2C_vertical(:, n_comb_v) = append_v;
                end
            end
            if y1_v >= fprop_sat_min - eps_v && y1_v <= fprop_sat_max + eps_v
                if z01_v >= fprop_sat_min - eps_v && z01_v <= fprop_sat_max + eps_v
                    append_v = [x0_v; y1_v];
                    n_comb_v = n_comb_v + 1;
                    P2C_vertical(:, n_comb_v) = append_v;
                end
            end
            if i_v > n_comb_v
                compute_v = false;
            end
        else
            if i_v < n_comb_v
                i_v = i_v+1;
            else
                compute_v = false;
            end
        end
    end
    P2C_vertical = P2C_vertical(:, 1:n_comb_v);
    n_v = size(P2C_vertical, 2);

    F_trans_computed_v = zeros(6,n_v);
    F3_test_v_w = zeros([1,n_v]);

    for i_v = 1:n_v
        F_trans_computed_v(:, i_v) = AutoFunc_sym_Fb_trans_w2(alpha_value, yaw, P2C_vertical(1,i_v), P2C_vertical(2,i_v))*~compute_via_w3 + AutoFunc_sym_Fb_trans_w3(alpha_value, yaw, P2C_vertical(1,i_v), P2C_vertical(2,i_v))*compute_via_w3;
        F3_test_v_w(i_v) = AutoFunc_sym_F3_trans_w2(alpha_value, yaw, P2C_vertical(1,i_v), P2C_vertical(2,i_v))*~compute_via_w3 + AutoFunc_sym_F3_trans_w3(alpha_value, yaw, P2C_vertical(1,i_v), P2C_vertical(2,i_v))*compute_via_w3;
    end

    if yaw == pi/2 || yaw == -pi/2
        Ft_v_num = F_trans_computed_v(2,:);
    else
        Ft_v_num = F_trans_computed_v(1,:)/cos(yaw);
    end

    Fz_v_num = F_trans_computed_v(3,:);

    [tri_verticale, ~] = convhull(Ft_v_num, Fz_v_num, "Simplify", true);

    x_cone_coupe = [0, r_cone, 0, -r_cone, 0];
    y_cone_coupe = [Fz_min, Fz_med, Fz_max, Fz_med, Fz_min];

    [~, k_bot] = min(Fz_v_num);
    [~, k_top] = max(Fz_v_num);
    [~, k_mid] = max(Ft_v_num);

    a_side_top = (Fz_v_num(k_mid)-Fz_v_num(k_top))/(Ft_v_num(k_mid)-Ft_v_num(k_top));
    a_side_bot = (Fz_v_num(k_mid)-Fz_v_num(k_bot))/(Ft_v_num(k_mid)-Ft_v_num(k_bot));
    b_side_top = Fz_v_num(k_top);
    b_side_bot = Fz_v_num(k_bot);

    side_top = @(t) a_side_top*t + b_side_top;
    side_bot = @(t) a_side_bot*t + b_side_bot;
    sphere_t = @(t) real((r_sphere.^2-t.^2).^(1/2));

    flag_Fd =  side_top(Fd_t)^2+Fd_t^2 >= r_sphere^2 && side_bot(Fd_t)^2 + Fd_t^2 <= r_sphere^2;

    if flag_Fd
        Fc_t = Fd_t;
        Fc_z = Fd(3);
        Fc = Fd;
        pitch_c = 0;
        % ratio = 1;
    else
        if alpha_value ~= 0 && r_sphere <= Fz_max && r_sphere >= Fz_min
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
            y0 = sphere_t(x0);
        else
            y0 = max([min([Fz_max, r_sphere]), Fz_min]);
        end
    
        %% Visus Rotations possibles Fd
    
        Fc_t = x0;
        Fc_z = y0;
        Fc = [x0*cos(yaw); x0*sin(yaw); y0];
    
        pitch_c = acos((Fd.'* Fc)/(norm(Fd)*norm(Fc)));
        % ratio = norm(Fc)/norm(Fd);
    end
else
    x0 = 0;
    y0 = max([min([Fz_max, r_sphere]), Fz_min]);
    Fc_t = x0;
    Fc_z = y0;
    Fc = [x0*cos(yaw); x0*sin(yaw); y0];
    pitch_c = acos((Fd.'* Fc)/(norm(Fd)*norm(Fc)));
    % ratio = norm(Fc)/norm(Fd);
end

%% PLOTS

Fx_limits = Plot_Limits(1,:);
Fy_limits = Plot_Limits(2,:);
Fz_limits = Plot_Limits(3,:);
Ft_limits = Plot_Limits(4,:);

frame = figure();
frame.Visible = 'off';
tiledlayout(1, 2);

nexttile();
hold on
if alpha_value~=0
    trisurf(tri_enveloppe, F_x_enveloppe, F_y_enveloppe, F_z_enveloppe, 'FaceColor', 'c', 'LineStyle', "-", 'EdgeColor', 'k', 'FaceAlpha', 0.2);
    surf(x_cone, y_cone, z_cone, 'FaceColor', 'y', 'LineStyle', 'none', 'FaceAlpha', 0.4);
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
    legend("Enveloppe forces atteignables", "Double cone simplifié des forces atteignables", ...
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
    plot(x_cone_coupe, y_cone_coupe, "y");
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
    legend("Forces atteignables pour $\psi="+num2str(yaw*180/pi)+"$°", "Coupe du cone simplifié", "$F_d$", "Cercle des rotations possibles de $F_d$", ...
        "Rotation optimale de $F_d$ de "+num2str(pitch_c*180/pi)+"°", 'Interpreter', 'latex', 'Location', 'northeast')
else
    legend("$F_d$", "Cercle des rotations possibles de $F_d$", ...
        "Rotation optimale de $F_d$ de "+num2str(pitch_c*180/pi)+"°", 'Interpreter', 'latex', 'Location', 'northeast')
end
xlim(Ft_limits);
ylim(Fz_limits);


end