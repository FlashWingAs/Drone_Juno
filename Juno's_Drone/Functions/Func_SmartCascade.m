function [Fc, Rot_Mat, ratio] = Func_SmartCascade(Fd, Mixer, alpha, sat_vec, dead_zone_alpha, alphaMesh_memory, yawMesh_memory, kFmot_memory)

eps_precision = 1e-5;
n_actionneur = 8;
fprop_sat_min = sat_vec(1) + 0.0*sat_vec(2);
fprop_sat_max = sat_vec(2) - 0.0*sat_vec(2);

Fz_max = (Mixer*fprop_sat_max*ones([n_actionneur,1])).'*[0; 0; 1; 0; 0; 0];
Fz_min = (Mixer*fprop_sat_min*ones([n_actionneur,1])).'*[0; 0; 1; 0; 0; 0];
% Fz_med = (Fz_max+Fz_min)/2;
% Try_vec = linspace(fprop_sat_min, fprop_sat_max, 2);

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

if alpha < dead_zone_alpha
    alpha_value = 0;
    cascade = true;
elseif alpha > pi - dead_zone_alpha
    alpha_value = pi;
    cascade = true;
else
    alpha_value = alpha;
    cascade = false;
end


%% Calcul slice verticale plan Fd

% Comb = Func_combinaisons(Try_vec, 2);
r_sphere2 = Fd_t^2 + Fd(3)^2;
r_sphere = sqrt(r_sphere2);

if ~cascade

    % [~, ~, kFmot, ~] = Func_SmartCascade_PreCompute(alpha_value, yaw, n_actionneur);
    % TO BE REPLACED WITH INTERPOLATION
    kFmot = Func_interpolationPreComputed(alpha_value, yaw, alphaMesh_memory, yawMesh_memory, kFmot_memory, dead_zone_alpha);

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

    % szFt = 0;
    % for i = 1:size(listFt, 2)
    %     if listFt(i)>=-eps_precision
    %         szFt = szFt + 1;
    %     end
    % end
    % listFt_compute = zeros(1,3);
    % listFz_compute = zeros(1,3);
    % 
    % k = 1;
    % for i = 1:size(listFt, 2)
    %     if listFt(i)>=-eps_precision
    %         listFt_compute(k) = listFt(i);
    %         listFz_compute(k) = listFz(i);
    %         k = k+1;
    %     end
    % end

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

Rot_Mat = Func_rot_z(yaw)*Func_rot_y(pitch_c)*Func_rot_z(-yaw);

end