
%% Creating interpolation grid

resolutionYaw   = pi/160;

interpolationAlpha = linspace(dead_zone_alpha, pi-dead_zone_alpha, 4);
interpolationYaw   = linspace(-pi+pi/720, pi, 720);

sizeAlpha = size(interpolationAlpha, 2);
sizeYaw   = size(interpolationYaw, 2);

% [alphaMesh_memory, yawMesh_memory] = meshgrid(interpolationAlpha, interpolationYaw);

flag_recompute = false;
if isfile("kFmot_memorySave.mat")
    load("kFmot_memorySave.mat", "kFmot_memory");
    if size(kFmot_memory,1) ~= sizeYaw || size(kFmot_memory,2) ~= sizeAlpha
        flag_recompute = true;
    end

else
    flag_recompute = true;
end

if flag_recompute
    kFmot_memory = zeros(sizeYaw, sizeAlpha, 2, 6);
    numericalErrorTreshold = 1e10;

    %% Computing every point on the grid
    % current_step_precompute = 0;
    % number_of_steps_precompute = sizeYaw*sizeAlpha;
    % objWaitBar_precompute = waitbar(current_step_precompute/number_of_steps_precompute, "pré-calcul en cours du controlleur - "+num2str(current_step_precompute/number_of_steps_precompute*100,3)+"%");

    parfor iAlpha = 1:sizeAlpha
        for iYaw = 1:sizeYaw
            temp = Func_SmartCascade_PreComputeLight(interpolationAlpha(iAlpha), interpolationYaw(iYaw));
            for i = 1:2
                for j = 1:6
                    if temp(i,j) > numericalErrorTreshold ||...
                            temp(i,j) < -numericalErrorTreshold
                        kFmot_memory(iYaw, iAlpha,i, j) = 0;
                    else
                        kFmot_memory(iYaw, iAlpha,i, j) = temp(i,j)
                    end
                end
            end
            % current_step_precompute = current_step_precompute + 1;
            % waitbar(current_step_precompute/number_of_steps_precompute, objWaitBar_precompute, "pré-calcul en cours du controlleur - "+num2str(current_step_precompute/number_of_steps_precompute*100,3)+"%")
        end
    end

    save('kFmot_memorySave.mat', 'kFmot_memory')
    close(objWaitBar_precompute)

end