Fd = [0; 0; m_drone_num*g_ref];
alpha_value = 30*pi/180;
dead_zone_alpha = 1*pi/180;
Plot_Limits = [-7.5, 7.5; -7.5, 7.5; -0, 30; -7.5, 7.5];

eps_precision = 1e-5;
n_actionneur = 8;
fprop_sat_min = fprop_sat_vec(1) + 0.1*fprop_sat_vec(2);
fprop_sat_max = fprop_sat_vec(2) - 0.0*fprop_sat_vec(2);
N_precision_sphere_cercle = 50;
t = linspace(0, 2*pi, N_precision_sphere_cercle);
if alpha_value < dead_zone_alpha
    alpha_value = 0;
end

Mixer = Func_sym_Mixer_wrapper(alpha_value);

Fz_max = (Mixer*fprop_sat_max*ones([n_actionneur,1])).'*[0; 0; 1; 0; 0; 0];
Fz_min = (Mixer*fprop_sat_min*ones([n_actionneur,1])).'*[0; 0; 1; 0; 0; 0];
Fz_med = (Fz_max+Fz_min)/2;
Try_vec = linspace(fprop_sat_min, fprop_sat_max, 2);

if Fd(1)~=0 || Fd(2)~=0
    yaw = atan2(Fd(2), Fd(1));
else
    yaw = 0;
end

if yaw == pi/2
    Fd_t = Fd(2);
elseif yaw == -pi/2
    Fd_t = -Fd(2);
else
    Fd_t = Fd(1)/cos(yaw);
end

if alpha_value < dead_zone_alpha
    alpha_value = 0;
    cascade = true;
elseif alpha_value > pi/2 - dead_zone_alpha
    alpha_value = pi;
    cascade = true;
else
    cascade = false;
end

%% Calcul enveloppe

if ~cascade

    P2C_enveloppe = Func_combinaisons(Try_vec, n_actionneur);
    n_enveloppe = size(P2C_enveloppe, 2);
    Enveloppe_drone = zeros(6, n_enveloppe);

    for i_e=1:n_enveloppe
        Enveloppe_drone(:, i_e) = Mixer*[P2C_enveloppe(:, i_e)];
    end

    F_x_enveloppe = Enveloppe_drone(1, :);
    F_y_enveloppe = Enveloppe_drone(2, :);
    F_z_enveloppe = Enveloppe_drone(3, :);

    [tri_enveloppe, ~] = convhull(F_x_enveloppe, F_y_enveloppe, F_z_enveloppe, "Simplify", false);

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

Comb = Func_combinaisons(Try_vec, 2);
r_sphere2 = Fd_t^2 + Fd(3)^2;
r_sphere = sqrt(r_sphere2);

if ~cascade

    % [~, ~, kFmot, ~] = Func_SmartCascade_PreCompute(alpha_value, yaw, n_actionneur);
    % TO BE REPLACED WITH INTERPOLATION
    kFmot = Func_interpolationPreComputed(alpha_value, yaw, interpolationAlpha, interpolationYaw, kFmot_memory, dead_zone_alpha);

    % Calcul des droites délimitants l'espace accessible pour F_mot

    % [a1; a2] * t + [b1; b2]
    % edgeLines{1, :} : cas pour sat-
    % edgeLines{2, :} : cas pour sat+
    edgeLines = zeros([2,n_actionneur,2,2]); 
    for i = 1:n_actionneur
        if ~Func_equal_eps(kFmot(2, i), 0, eps_precision)
            edgeLines(1, i, :, :) = [1, 0; -kFmot(1,i)/kFmot(2,i), fprop_sat_min/kFmot(2,i)];
            edgeLines(2, i, :, :) = [1, 0; -kFmot(1,i)/kFmot(2,i), fprop_sat_max/kFmot(2,i)];
        elseif ~Func_equal_eps(kFmot(1, i), 0, eps_precision)
            edgeLines(1, i, :, :) = [-kFmot(2,i)/kFmot(1,i), fprop_sat_min/kFmot(1,i); 1, 0];
            edgeLines(2, i, :, :) = [-kFmot(2,i)/kFmot(1,i), fprop_sat_max/kFmot(1,i); 1, 0];
        end
    end

    % Calcul des intersections des droites
    % nombre max d'intersections = 2n C 2
    listIntersectionsMax = Inf([2, (2*n_actionneur)*(2*n_actionneur-1)/2]);
    k = 1;
    for i = 1 : n_actionneur-1
        for j = i+1 : n_actionneur
            listIntersectionsMax(:,k+0) = Func_intersect_parametrized(edgeLines(1,j,:,1), edgeLines(1,j,:,2), edgeLines(1,i,:,1), edgeLines(1,i,:,2));
            listIntersectionsMax(:,k+1) = Func_intersect_parametrized(edgeLines(2,j,:,1), edgeLines(2,j,:,2), edgeLines(2,i,:,1), edgeLines(2,i,:,2));
            listIntersectionsMax(:,k+2) = Func_intersect_parametrized(edgeLines(1,j,:,1), edgeLines(1,j,:,2), edgeLines(2,i,:,1), edgeLines(2,i,:,2));
            listIntersectionsMax(:,k+3) = Func_intersect_parametrized(edgeLines(2,j,:,1), edgeLines(2,j,:,2), edgeLines(1,i,:,1), edgeLines(1,i,:,2));
            k = k + 4;
        end
    end
    
    [~, listIntersections_col] = find(listIntersectionsMax ~= Inf);
    listIntersections_raw =  Func_remDuplicate(listIntersectionsMax(:, listIntersections_col));
    indices2keep = true([1, size(listIntersections_raw, 2)]);
    for i = 1:size(listIntersections_raw, 2)
        for k = 1:n_actionneur
            temp = kFmot(1, k)*listIntersections_raw(1, i) + kFmot(2, k)*listIntersections_raw(2, i);
            if Func_greater_eps(temp, fprop_sat_max, eps_precision)
                indices2keep(i) = false;
                break;
            elseif Func_lower_eps(temp, fprop_sat_min, eps_precision)
                indices2keep(i) = false;
                break;
            end
        end
    end

    listIntersections = listIntersections_raw(:, indices2keep);
    nIntersection = size(listIntersections, 2);

    %% Passage des vertices des lambda aux vertices de Fmot
    listFmot = zeros([n_actionneur, nIntersection]);
    for i = 1:nIntersection
        for j = 1:n_actionneur
            listFmot(j, i) = kFmot(1, j)*listIntersections(1, i) + kFmot(2, j)*listIntersections(2, i);
        end
    end

    %% Passage aux vertices dans le plan (Ft, Fz)
    listFw = zeros([6, nIntersection]);
    for i = 1:nIntersection
        listFw(:,i) = Mixer*listFmot(:,i);        
    end

    if yaw == pi/2 || yaw == -pi/2
        listFt = listFw(2,:);
    else
        listFt = listFw(1,:)/cos(yaw);
    end
    listFz = listFw(3,:);

    %% Génération des droites dans les plan (Ft, Fz)

    [tri_plan_visu, ~] = convhull(listFt, listFz, "Simplify", false);
    listFt_visu = listFt(tri_plan_visu);
    listFz_visu = listFz(tri_plan_visu);
    nLines_visu = size(listFt_visu, 2)-1;

    %%% paramètres de la droite,
    %%% point de début, point de fin.
    listLinesFw_visu = cell([2, nLines_visu]);
    for i = 1:nLines_visu
        temp = zeros([2,2]);
        x1 = listFt_visu(i);
        y1 = listFz_visu(i);
        if i < nLines_visu
            x2 = listFt_visu(i+1);
            y2 = listFz_visu(i+1);
        else
            x2 = listFt_visu(1);
            y2 = listFz_visu(1);
        end
        if ~Func_equal_eps(x1, x2, eps_precision)
            temp(1,1) = 1;
            temp(1,2) = 0;
            temp(2,1) = (y2-y1)/(x2-x1);
            temp(2,2) = y1-temp(2,1)*x1;
        elseif ~Func_equal_eps(y1, y2, eps_precision)
            temp(1,1) = (x2-x1)/(y2-y1);
            temp(1,2) = x1-temp(1,1)*y1;
            temp(2,1) = 1;
            temp(2,2) = 0;
        %%% pas de else car 2 points ne peuvent être au même endroit, car
        %%% retrait des doublons avant
        end
        listLinesFw_visu{1,i} = temp;
        if i < nLines_visu
            listLinesFw_visu{2,i} = [Func_root_parmLine(temp(:,1), temp(:,2), [listFt_visu(i); listFz_visu(i)]), Func_root_parmLine(temp(:,1), temp(:,2), [listFt_visu(i+1); listFz_visu(i+1)])];
        else
            listLinesFw_visu{2,i} = [Func_root_parmLine(temp(:,1), temp(:,2), [listFt_visu(i); listFz_visu(i)]), Func_root_parmLine(temp(:,1), temp(:,2), [listFt_visu(1); listFz_visu(1)])];
        end
    end
    listFt_compute = listFt(listFt>=-eps_precision);
    listFz_compute = listFz(listFt>=-eps_precision);


    % [tri_plan_compute, ~] = convhull(listFt_compute, listFz_compute, "Simplify", false);
    [~, k_sort] = sort(listFz_compute);
    tri_plan_compute = [k_sort, k_sort(1)];
    
    
    listFt_compute = listFt_compute(tri_plan_compute);
    listFz_compute = listFz_compute(tri_plan_compute);
    nLines_compute = size(listFt_compute, 2)-2;

    %%% paramètres de la droite,
    %%% point de début, point de fin.
    %%% Fd à l'exterieur
    listLinesFw_computeLine = zeros([nLines_compute, 2, 2]);
    listLinesFw_computePoints = zeros([nLines_compute, 2]);
    listLinesFw_computeBool = false(nLines_compute);
    for i = 1:nLines_compute
        temp = zeros([2,2]);
        x1 = listFt_compute(i);
        y1 = listFz_compute(i);
        % if i < nLines_compute
            x2 = listFt_compute(i+1);
            y2 = listFz_compute(i+1);
        % else
            % x2 = listFt_visu(1);
            % y2 = listFz_visu(1);
        % end
        if ~Func_equal_eps(x1, x2, eps_precision)
            temp(1,1) = 1;
            temp(1,2) = 0;
            temp(2,1) = (y2-y1)/(x2-x1);
            temp(2,2) = y1-temp(2,1)*x1;
        elseif ~Func_equal_eps(y1, y2, eps_precision)
            temp(1,1) = (x2-x1)/(y2-y1);
            temp(1,2) = x1-temp(1,1)*y1;
            temp(2,1) = 1;
            temp(2,2) = 0;
        %%% pas de else car 2 points ne peuvent être au même endroit, car
        %%% retrait des doublons avant
        end
        listLinesFw_computeLine(i,:,:) = temp;
        % if i < nLines_compute
            listLinesFw_computePoints(i,:) = [Func_root_parmLine(temp(:,1), temp(:,2), [listFt_compute(i); listFz_compute(i)]), Func_root_parmLine(temp(:,1), temp(:,2), [listFt_compute(i+1); listFz_compute(i+1)])];
        % else
            % listLinesFw_compute{2,i} = [Func_root_parmLine(temp(:,1), temp(:,2), [listFt_compute(i); listFz_compute(i)]), Func_root_parmLine(temp(:,1), temp(:,2), [listFt_compute(1); listFz_compute(1)])];
        % end
        listLinesFw_computeBool(i) = Func_RightOrLeft([x2-x1; y2-y1], [Fd_t-x1; Fd(3)-y1]);
    end

    %% Vérification si Fd est en dehors de l'espace de forces accessibles

    flag_Fd = true;
    for i = 1:nLines_compute
        if listLinesFw_computeBool(i)
            flag_Fd = false;
        end
    end
    
    sphere_t = @(t) real((r_sphere.^2-t.^2).^(1/2));

    if flag_Fd
        Fc_t = Fd_t;
        Fc_z = Fd(3);
        Fc = Fd;
        pitch_c = 0;
        ratio = 1;
    else
        if alpha_value ~= 0 && r_sphere <= Fz_max && r_sphere >= Fz_min
            x0 = Inf;
            for i = 1:nLines_compute
                poly = [listLinesFw_computeLine(i,1,1).^2 + listLinesFw_computeLine(i,2,1).^2, ...
                        2*listLinesFw_computeLine(i,1,1)*listLinesFw_computeLine(i,1,2) + 2*listLinesFw_computeLine(i,2,1)*listLinesFw_computeLine(i,2,2), ...
                        listLinesFw_computeLine(i,1,2).^2 + listLinesFw_computeLine(i,2,2).^2 - r_sphere2];
                roots_temp = roots(poly);
                x_temp = min(abs(roots_temp));
                if x0 == Inf
                    x0 = x_temp;
                elseif x_temp < x0
                    x0 = x_temp;
                end                
            end

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
        ratio = norm(Fc)/norm(Fd);
    end
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
ticks_horizontal = [-7.5, -5, -2.5, 0, 2.5, 5, 7.5];

fontSize_title = 30;
fontSize_legend = 15;
fontSize_axisName = 30;
fontSize_axisTicks = 20;

frame = figure();
frame.Visible = 'off';
tiledlayout(1, 3);

tile = nexttile();
hold on
if ~cascade
    trisurf(tri_enveloppe, F_x_enveloppe, F_y_enveloppe, F_z_enveloppe, 'FaceColor', 'c', 'LineStyle', "-", 'LineWidth', 1.5, 'EdgeColor', 'k', 'FaceAlpha', 0.5, 'DisplayName', "Enveloppe forces atteignables");
end
plot3([0, Fd(1)], [0, Fd(2)], [0, Fd(3)], '-rx', 'LineWidth', 2, 'MarkerSize', 20, 'DisplayName', "$F_d$");
surf(x_sphere, y_sphere, z_sphere, 'FaceColor', 'r', 'LineStyle', 'none', 'FaceAlpha', 0.5, 'DisplayName', "Sphere des rotations possibles de $F_d$")
plot3([0, Fc(1)], [0, Fc(2)], [0, Fc(3)], '--go', 'LineWidth', 2, 'MarkerSize', 20, 'DisplayName', "Rotation optimale de $F_d$ de "+num2str(pitch_c*180/pi)+"°")
hold off
axis equal
grid();
xlabel('$F_x$', 'Interpreter', 'latex', 'FontSize', fontSize_axisName);
ylabel('$F_y$', 'Interpreter', 'latex', 'FontSize', fontSize_axisName);
zlabel('$F_z$', 'Interpreter', 'latex', 'FontSize', fontSize_axisName, 'Rotation', 0);
title("Volume des forces atteignables"+newline+"pour $\alpha=$"+num2str(alpha_value*180/pi)+"$^\circ$", 'Interpreter', 'latex', 'FontSize', fontSize_title)
legend('Interpreter', 'latex', 'Location', 'northeast', 'FontSize', fontSize_legend)
view(45, 30)
lightangle(-45,30)
lighting flat
xlim(Fx_limits);
ylim(Fy_limits);
zlim(Fz_limits);
xticks(ticks_horizontal)
yticks(ticks_horizontal)
xtickangle(0)
ytickangle(0)
tile.XAxis.FontSize = fontSize_axisTicks;
tile.YAxis.FontSize = fontSize_axisTicks;
tile.ZAxis.FontSize = fontSize_axisTicks;

tile = nexttile();
hold on
if ~cascade
    t = -100:0.1:100;
    for i = 1:n_actionneur
        d1 = reshape(edgeLines(1,i,:,1), [2,1])*t+reshape(edgeLines(1,i,:,2), [2,1]);
        d2 = reshape(edgeLines(2,i,:,1), [2,1])*t+reshape(edgeLines(2,i,:,2), [2,1]);
        plot(d1(1,:), d1(2,:), '-b', 'LineWidth', 2, 'HandleVisibility','off')
        plot(d2(1,:), d2(2,:), '-r', 'LineWidth', 2, 'HandleVisibility','off')
    end
    % for i = 1:size(listIntersections, 2)
    scatter(listIntersections_raw(1,:), listIntersections_raw(2,:), 100, 'white', 'o', "filled", 'DisplayName', "Edges intersections")
    % end
    [tri_scatter, ~] = convhull(listIntersections(1,:), listIntersections(2,:));
    fill(listIntersections(1,tri_scatter), listIntersections(2,tri_scatter), 'cyan', 'FaceColor', 'cyan', 'EdgeColor', 'cyan', 'FaceAlpha', 0.5, 'DisplayName', "Espace born\'e par les limites")

end
xlim([fprop_sat_min - 0.5*(fprop_sat_max-fprop_sat_min), fprop_sat_max + 0.5*(fprop_sat_max-fprop_sat_min)])
ylim([fprop_sat_min - 0.5*(fprop_sat_max-fprop_sat_min), fprop_sat_max + 0.5*(fprop_sat_max-fprop_sat_min)])
grid()
xlabel('$\lambda_1$', 'Interpreter', 'latex', 'FontSize', fontSize_axisName);
ylabel('$\lambda_2$', 'Interpreter', 'latex', 'FontSize', fontSize_axisName, 'Rotation', 0);
title("Plan ($\lambda_1, \lambda_2$),"+newline+"d\'elimit\'e par les saturations", 'Interpreter', 'latex', 'FontSize', fontSize_title)
legend('Interpreter', 'latex', 'Location', 'northeast', 'FontSize', fontSize_legend)
tile.XAxis.FontSize = fontSize_axisTicks;
tile.YAxis.FontSize = fontSize_axisTicks;
tile.ZAxis.FontSize = fontSize_axisTicks;


tile = nexttile();
hold on
if ~cascade
    for i=1:nLines_visu
        t = listLinesFw_visu{2,i};
        d = listLinesFw_visu{1,i}(:,1)*t+listLinesFw_visu{1,i}(:,2);
        if i>1
            plot(d(1,:), d(2,:), "c", 'LineWidth', 2, 'HandleVisibility','off');
        else
            plot(d(1,:), d(2,:), "c", 'LineWidth', 2, 'DisplayName', "Forces atteignables pour $\psi="+num2str(yaw*180/pi)+"$°");
        end
    end
end
plot([0, Fd_t], [0, Fd(3)], '-rx', 'LineWidth', 2, 'MarkerSize', 20, 'DisplayName', "$F_d$")
plot(x_circle_force, y_circle_force, '--r', 'LineWidth', 2, 'DisplayName', "Cercle des rotations possibles de $F_d$")
plot([0, Fc_t], [0, Fc_z], '--go', 'LineWidth', 2, 'MarkerSize', 20, 'DisplayName', "Rotation optimale de $F_d$ de "+num2str(pitch_c*180/pi)+"°")
hold off
axis equal
grid()
xlabel('$F_t$', 'Interpreter', 'latex', 'FontSize', fontSize_axisName);
ylabel('$F_z$', 'Interpreter', 'latex', 'FontSize', fontSize_axisName, 'Rotation', 0);
title("Volume des forces atteignables pour $\alpha=$"+num2str(alpha_value*180/pi)+"$^\circ$,"+newline+"dans le plan ($F_t,F_z)$", 'Interpreter', 'latex', 'FontSize', fontSize_title)
legend('Interpreter', 'latex', 'Location', 'northeast', 'FontSize', fontSize_legend)
xlim(Ft_limits);
ylim(Fz_limits);
xticks(ticks_horizontal)
xtickangle(0)
tile.XAxis.FontSize = fontSize_axisTicks;
tile.YAxis.FontSize = fontSize_axisTicks;
tile.ZAxis.FontSize = fontSize_axisTicks;


frame.Visible = 'on';