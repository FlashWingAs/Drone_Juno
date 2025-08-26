%% Données Simulink

text = 'Récupération des données';
Func_Progress_Log(false, text, LOG_SIZE);

% loggedData = simOut1.logsout;
loggedData = out.logsout;

Func_Progress_Log(true, text, LOG_SIZE);

%% Traitement des données

text = 'Traitement des données';
Func_Progress_Log(false, text, LOG_SIZE);

ts_p_mes = loggedData.get("p0_mes").Values;
ts_dp_mes = loggedData.get("dp0_mes").Values;
ts_ddp_mes = loggedData.get("ddp0_mes").Values;
ts_eta_mes = loggedData.get("eta_mes").Values;
ts_R0B_mes = loggedData.get("R0B_mes").Values;
ts_omega_mes = loggedData.get("omega_mes").Values;
ts_domega_mes = loggedData.get("domega_mes").Values;

p_mes = reshape(ts_p_mes.Data(:,:,:), [3, numel(ts_p_mes.Data(1,:,:))]);
px_mes = reshape(ts_p_mes.Data(1,:,:), [numel(ts_p_mes.Data(1,:,:)), 1]);
py_mes = reshape(ts_p_mes.Data(2,:,:), [numel(ts_p_mes.Data(2,:,:)), 1]);
pz_mes = reshape(ts_p_mes.Data(3,:,:), [numel(ts_p_mes.Data(3,:,:)), 1]);

dp_mes = reshape(ts_dp_mes.Data(:,:,:), [3, numel(ts_dp_mes.Data(1,:,:))]);
dpx_mes = reshape(ts_dp_mes.Data(1,:,:), [numel(ts_dp_mes.Data(1,:,:)), 1]);
dpy_mes = reshape(ts_dp_mes.Data(2,:,:), [numel(ts_dp_mes.Data(2,:,:)), 1]);
dpz_mes = reshape(ts_dp_mes.Data(3,:,:), [numel(ts_dp_mes.Data(3,:,:)), 1]);

ddpx_mes = reshape(ts_ddp_mes.Data(1,:,:), [numel(ts_ddp_mes.Data(1,:,:)), 1]);
ddpy_mes = reshape(ts_ddp_mes.Data(2,:,:), [numel(ts_ddp_mes.Data(2,:,:)), 1]);
ddpz_mes = reshape(ts_ddp_mes.Data(3,:,:), [numel(ts_ddp_mes.Data(3,:,:)), 1]);

etax_mes = reshape(ts_eta_mes.Data(1,:,:), [numel(ts_eta_mes.Data(1,:,:)), 1]);
etay_mes = reshape(ts_eta_mes.Data(2,:,:), [numel(ts_eta_mes.Data(2,:,:)), 1]);
etaz_mes = reshape(ts_eta_mes.Data(3,:,:), [numel(ts_eta_mes.Data(3,:,:)), 1]);

R0B_mes = reshape(ts_R0B_mes.Data(:,:,:), [3, 3, numel(ts_R0B_mes.Data(3,3,:))]);

omegax_mes = reshape(ts_omega_mes.Data(1,:,:), [numel(ts_omega_mes.Data(1,:,:)), 1]);
omegay_mes = reshape(ts_omega_mes.Data(2,:,:), [numel(ts_omega_mes.Data(2,:,:)), 1]);
omegaz_mes = reshape(ts_omega_mes.Data(3,:,:), [numel(ts_omega_mes.Data(3,:,:)), 1]);

domegax_mes = reshape(ts_domega_mes.Data(1,:,:), [numel(ts_domega_mes.Data(1,:,:)), 1]);
domegay_mes = reshape(ts_domega_mes.Data(2,:,:), [numel(ts_domega_mes.Data(2,:,:)), 1]);
domegaz_mes = reshape(ts_domega_mes.Data(3,:,:), [numel(ts_domega_mes.Data(3,:,:)), 1]);

ts_p_com = loggedData.get("<p0_com>").Values;
ts_dp_com = loggedData.get("<dp0_com>").Values;
ts_ddp_com = loggedData.get("<ddp0_com>").Values;
ts_eta_com = loggedData.get("<eta_com>").Values;
ts_R0B_com = loggedData.get("<R0B_com>").Values;
ts_omega_com = loggedData.get("<omega_com>").Values;
ts_domega_com = loggedData.get("<domega_com>").Values;
ts_alpha_com = loggedData.get("alpha_com").Values;

p_com = reshape(ts_p_com.Data(:,:,:), [3, numel(ts_p_com.Data(1,:,:))]);
px_com = reshape(ts_p_com.Data(1,:,:), [numel(ts_p_com.Data(1,:,:)), 1]);
py_com = reshape(ts_p_com.Data(2,:,:), [numel(ts_p_com.Data(2,:,:)), 1]);
pz_com = reshape(ts_p_com.Data(3,:,:), [numel(ts_p_com.Data(3,:,:)), 1]);

dp_com = reshape(ts_dp_com.Data(:,:,:), [3, numel(ts_dp_com.Data(1,:,:))]);
dpx_com = reshape(ts_dp_com.Data(1,:,:), [numel(ts_dp_com.Data(1,:,:)), 1]);
dpy_com = reshape(ts_dp_com.Data(2,:,:), [numel(ts_dp_com.Data(2,:,:)), 1]);
dpz_com = reshape(ts_dp_com.Data(3,:,:), [numel(ts_dp_com.Data(3,:,:)), 1]);

ddp_com = reshape(ts_ddp_com.Data(:,:,:), [3, numel(ts_ddp_com.Data(1,:,:))]);
ddpx_com = reshape(ts_ddp_com.Data(1,:,:), [numel(ts_ddp_com.Data(1,:,:)), 1]);
ddpy_com = reshape(ts_ddp_com.Data(2,:,:), [numel(ts_ddp_com.Data(2,:,:)), 1]);
ddpz_com = reshape(ts_ddp_com.Data(3,:,:), [numel(ts_ddp_com.Data(3,:,:)), 1]);

etax_com = reshape(ts_eta_com.Data(1,:,:), [numel(ts_eta_com.Data(1,:,:)), 1]);
etay_com = reshape(ts_eta_com.Data(2,:,:), [numel(ts_eta_com.Data(2,:,:)), 1]);
etaz_com = reshape(ts_eta_com.Data(3,:,:), [numel(ts_eta_com.Data(3,:,:)), 1]);

R0B_com = reshape(ts_R0B_com.Data(:,:,:), [3, 3, numel(ts_R0B_com.Data(3,3,:))]);

omegax_com = reshape(ts_omega_com.Data(1,:,:), [numel(ts_omega_com.Data(1,:,:)), 1]);
omegay_com = reshape(ts_omega_com.Data(2,:,:), [numel(ts_omega_com.Data(2,:,:)), 1]);
omegaz_com = reshape(ts_omega_com.Data(3,:,:), [numel(ts_omega_com.Data(3,:,:)), 1]);

domegax_com = reshape(ts_domega_com.Data(1,:,:), [numel(ts_domega_com.Data(1,:,:)), 1]);
domegay_com = reshape(ts_domega_com.Data(2,:,:), [numel(ts_domega_com.Data(2,:,:)), 1]);
domegaz_com = reshape(ts_domega_com.Data(3,:,:), [numel(ts_domega_com.Data(3,:,:)), 1]);

alpha_com = reshape(ts_alpha_com.Data(:,:), [numel(ts_alpha_com.Data(:,:)), 1]);

ts_Time = ts_p_mes.Time;
n_ts = numel(ts_Time);

e_x = p_com-p_mes;
e_v = dp_com-dp_mes;
Fr_base = zeros([3, n_ts]);

for i = 1:n_ts
    Fr_base(:,i) = R0B_com(:,:,i).'*(k_corr_x*e_x(:,i) + k_corr_v*e_v(:,i) + m_drone_num*g_ref*[0;0;1] + m_drone_num*ddp_com(:,i));
end


Func_Progress_Log(true, text, LOG_SIZE);

%% Film des efforts

film = true;
if film
    init_frame_rate = 1000;
    desired_frame_rate = 30;
    writerObj = VideoWriter('Video_efforts_accessibles_Cercle.mp4', 'MPEG-4');
    writerObj.FrameRate = desired_frame_rate;
    open(writerObj)
    Lims = [-10, 10; -10, 10; -2, 25; -10, 10];
    number_of_frames = floor(n_ts*desired_frame_rate/init_frame_rate);
    frames = cell(1,number_of_frames);
    f = waitbar(0, "0% --- 0 images sur "+num2str(number_of_frames)+" traitées");
    for i = 1:number_of_frames
        frames{i} = getframe(Func_visu_movie(Fr_base(:,ceil(i*init_frame_rate/desired_frame_rate)), alpha_com(ceil(i*init_frame_rate/desired_frame_rate)), dead_zone_alpha, fprop_sat_vec, Lims));
        
        f = waitbar(i/number_of_frames, f, num2str(i/number_of_frames*100,5)+" % --- "+num2str(i)+" images sur "+num2str(number_of_frames)+" traitées");
    end
    close(f);
    disp("frames crées")
    f = waitbar(0, "0% --- génération vidéo");
    for i = 1:number_of_frames
        writeVideo(writerObj, frames{i});
        f = waitbar(i/number_of_frames, f, "génération vidéo");
    end
    close(writerObj);
    close(f);
    disp("vidéo crée")
end


%% Comparaisons directes

if ~film
    f1 = figure();
    t1 = tiledlayout(3,4);

    nexttile()
    grid();
    hold on
    plot(ts_Time, px_com);
    plot(ts_Time, px_mes);
    fill([ts_Time;flipud(ts_Time)], [px_com;flipud(px_mes)], [0.9 0.9 1], 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    ylabel('$p_x(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    legend('Commande', 'Mesures', 'Location','northeast');
    hold off

    nexttile()
    grid();
    hold on
    plot(ts_Time, dpx_com);
    plot(ts_Time, dpx_mes);
    fill([ts_Time;flipud(ts_Time)], [dpx_com;flipud(dpx_mes)], [0.9 0.9 1], 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    ylabel('$\dot{p}_x(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    legend('Commande', 'Mesures', 'Location','northeast');
    hold off

    nexttile()
    grid();
    hold on
    plot(ts_Time, etax_com);
    plot(ts_Time, etax_mes);
    fill([ts_Time;flipud(ts_Time)], [etax_com;flipud(etax_mes)], [0.9 0.9 1], 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    ylabel('$\eta_x(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    legend('Commande', 'Mesures', 'Location','northeast');
    hold off

    nexttile()
    grid();
    hold on
    plot(ts_Time, omegax_com);
    plot(ts_Time, omegax_mes);
    fill([ts_Time;flipud(ts_Time)], [omegax_com;flipud(omegax_mes)], [0.9 0.9 1], 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    ylabel('$\omega_x(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    legend('Commande', 'Mesures', 'Location','northeast');
    hold off

    nexttile()
    grid();
    hold on
    plot(ts_Time, py_com);
    plot(ts_Time, py_mes);
    fill([ts_Time;flipud(ts_Time)], [py_com;flipud(py_mes)], [0.9 0.9 1], 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    ylabel('$p_y(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    legend('Commande', 'Mesures', 'Location','northeast');
    hold off

    nexttile()
    grid();
    hold on
    plot(ts_Time, dpy_com);
    plot(ts_Time, dpy_mes);
    fill([ts_Time;flipud(ts_Time)], [dpy_com;flipud(dpy_mes)], [0.9 0.9 1], 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    ylabel('$\dot{p}_y(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    legend('Commande', 'Mesures', 'Location','northeast');
    hold off

    nexttile()
    grid();
    hold on
    plot(ts_Time, etay_com);
    plot(ts_Time, etay_mes);
    fill([ts_Time;flipud(ts_Time)], [etay_com;flipud(etay_mes)], [0.9 0.9 1], 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    ylabel('$\eta_y(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    legend('Commande', 'Mesures', 'Location','northeast');
    hold off

    nexttile()
    grid();
    hold on
    plot(ts_Time, omegay_com);
    plot(ts_Time, omegay_mes);
    fill([ts_Time;flipud(ts_Time)], [omegay_com;flipud(omegay_mes)], [0.9 0.9 1], 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    ylabel('$\omega_y(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    legend('Commande', 'Mesures', 'Location','northeast');
    hold off

    nexttile()
    grid();
    hold on
    plot(ts_Time, pz_com);
    plot(ts_Time, pz_mes);
    fill([ts_Time;flipud(ts_Time)], [pz_com;flipud(pz_mes)], [0.9 0.9 1], 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    ylabel('$p_z(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    legend('Commande', 'Mesures', 'Location','northeast');
    hold off

    nexttile()
    grid();
    hold on
    plot(ts_Time, dpz_com);
    plot(ts_Time, dpz_mes);
    fill([ts_Time;flipud(ts_Time)], [dpz_com;flipud(dpz_mes)], [0.9 0.9 1], 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    ylabel('$\dot{p}_z(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    legend('Commande', 'Mesures', 'Location','northeast');
    hold off

    nexttile()
    grid();
    hold on
    plot(ts_Time, etaz_com);
    plot(ts_Time, etaz_mes);
    fill([ts_Time;flipud(ts_Time)], [etaz_com;flipud(etaz_mes)], [0.9 0.9 1], 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    ylabel('$\eta_z(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    legend('Commande', 'Mesures', 'Location','northeast');
    hold off

    nexttile()
    grid();
    hold on
    plot(ts_Time, omegaz_com);
    plot(ts_Time, omegaz_mes);
    fill([ts_Time;flipud(ts_Time)], [omegaz_com;flipud(omegaz_mes)], [0.9 0.9 1], 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    ylabel('$\omega_z(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    legend('Commande', 'Mesures', 'Location','northeast');
    hold off

    sgtitle("Comparaisons directes")

    %% Ecarts

    f2 = figure();
    t2 = tiledlayout(3,4);

    nexttile()
    grid();
    hold on
    plot(ts_Time, px_com-px_mes);
    ylabel('$p_x(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    hold off

    nexttile()
    grid();
    hold on
    plot(ts_Time, dpx_com-dpx_mes);
    ylabel('$\dot{p}_x(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    hold off

    nexttile()
    grid();
    hold on
    plot(ts_Time, etax_com-etax_mes);
    ylabel('$\eta_x(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    hold off

    nexttile()
    grid();
    hold on
    plot(ts_Time, omegax_com-omegax_mes);
    ylabel('$\omega_x(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    hold off

    nexttile()
    grid();
    hold on
    plot(ts_Time, py_com-py_mes);
    ylabel('$p_y(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    hold off

    nexttile()
    grid();
    hold on
    plot(ts_Time, dpy_com-dpy_mes);
    ylabel('$\dot{p}_y(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    hold off

    nexttile()
    grid();
    hold on
    plot(ts_Time, etay_com-etay_mes);
    ylabel('$\eta_y(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    hold off

    nexttile()
    grid();
    hold on
    plot(ts_Time, omegay_com-omegay_mes);
    ylabel('$\omega_y(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    hold off

    nexttile()
    grid();
    hold on
    plot(ts_Time, pz_com-pz_mes);
    ylabel('$p_z(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    hold off

    nexttile()
    grid();
    hold on
    plot(ts_Time, dpz_com-dpz_mes);
    ylabel('$\dot{p}_z(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    hold off

    nexttile()
    grid();
    hold on
    plot(ts_Time, etaz_com-etaz_mes);
    ylabel('$\eta_z(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    hold off

    nexttile()
    grid();
    hold on
    plot(ts_Time, omegaz_com-omegaz_mes);
    ylabel('$\omega_z(t)$', 'Interpreter', 'latex', 'Rotation', 0);
    xlabel('Time(s)', 'Interpreter', 'latex');
    hold off

    sgtitle("Mesure des écarts")
end