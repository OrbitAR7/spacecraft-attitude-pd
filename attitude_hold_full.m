
%% ==========================================
%  SPACECRAFT ATTITUDE CONTROL SIMULATION 
%  File: main.m
%  - PD attitude hold with 4 reaction wheels
%  - Robust shapes, logging, metrics
%  - Generates multiple figures & saves PNGs
%% ==========================================

clear; clc; close all;

fprintf('==========================================\n');
fprintf('SPACECRAFT ATTITUDE CONTROL SIMULATION\n');
fprintf('==========================================\n\n');

%% ---- Paths (relative only; safe-guarded) ----
rel_paths = {'config','src/dynamics','src/control','src/utilities','src/visualization'};
for p = rel_paths
    if isfolder(p{1}), addpath(p{1}); end %#ok<AGROW>
end

set(0,'DefaultFigureVisible','on');  % ensure figures show

%% ---- Simulation parameters ----
sim_params.duration = 1000.0;    % seconds
sim_params.dt       = 0.1;       % seconds
sim_params.t        = 0:sim_params.dt:sim_params.duration;
n_steps = numel(sim_params.t);

%% ---- Spacecraft & wheels ----
spacecraft.inertia = diag([6400, 4730, 8160]);   % kg*m^2
spacecraft.quaternion = [1 0 0 0];               % [qw qx qy qz]
spacecraft.angular_velocity = [0.01 0.01 0.01];  % rad/s

waxes = [1 0 0;
         0 1 0;
         0 0 1;
         1 1 1];
waxes(4,:) = waxes(4,:)./norm(waxes(4,:));
wheel_inertia = 0.05;        % kg*m^2

spacecraft.num_wheels = size(waxes,1);
for i = 1:spacecraft.num_wheels
    spacecraft.wheels(i).axis    = waxes(i,:);  %#ok<SAGROW>
    spacecraft.wheels(i).inertia = wheel_inertia;
    spacecraft.wheels(i).speed   = 0.0;         % rad/s
end

%% ---- Controller ----
controller_params.Kp = 10.0;
controller_params.Kd = 450.0;
q_ref    = [1 0 0 0];
w_ref    = [0 0 0];

%% ---- Disturbance (set to zero) ----
disturbance = [0;0;0];

%% ---- History ----
history.t                = sim_params.t(:);
history.quaternion       = zeros(n_steps,4);
history.angular_velocity = zeros(n_steps,3);
history.wheel_speeds     = zeros(n_steps, spacecraft.num_wheels);
history.control_torque   = zeros(n_steps,3);
history.attitude_error   = zeros(n_steps,1);

history = store_state(history, spacecraft, 1);

fprintf('Simulation Configuration:\n');
fprintf('  Duration: %.1f seconds\n', sim_params.duration);
fprintf('  Time step: %.3f seconds\n', sim_params.dt);
fprintf('  Spacecraft inertia: [%g, %g, %g] kg*m^2\n', diag(spacecraft.inertia));
fprintf('  Number of reaction wheels: %d\n', spacecraft.num_wheels);
fprintf('  Controller gains: Kp=%.1f, Kd=%.1f\n\n', controller_params.Kp, controller_params.Kd);

fprintf('Starting simulation...\nProgress: \n');
fprintf('Initial state check:\n');
fprintf('  Quaternion size: %dx%d\n', size(spacecraft.quaternion,1), size(spacecraft.quaternion,2));
fprintf('  Angular velocity size: %dx%d\n', size(spacecraft.angular_velocity,1), size(spacecraft.angular_velocity,2));
fprintf('  Angular velocity values: [%.4f, %.4f, %.4f]\n\n', spacecraft.angular_velocity);

%% ---- Main loop ----
for k = 2:n_steps
    % Perfect estimation
    q_est = spacecraft.quaternion;
    w_est = spacecraft.angular_velocity;
    % Guidance & error
    error_state = compute_error_state(q_ref, q_est, w_ref, w_est);
    % Control
    control_torque = apply_control_law(error_state, controller_params);
    history.attitude_error(k)   = error_state.angle;
    history.control_torque(k,:) = control_torque;
    % Allocation
    wheel_torques = allocate_to_wheels(spacecraft, control_torque);
    % Dynamics
    spacecraft = propagate_dynamics(spacecraft, wheel_torques, disturbance, sim_params.dt);
    % Log
    history = store_state(history, spacecraft, k);
    % Progress
    if mod(k, max(1, round(0.05*n_steps))) == 0
        pct = round(100*k/n_steps);
        bars = repmat('=',1,round(pct/5));
        spaces = repmat(' ',1,20 - numel(bars));
        fprintf('Progress: [%s%s] %3d%%\n', bars, spaces, pct);
    end
end

fprintf('\nSimulation completed.\n\nPost-processing results...\n\n');

%% ---- Metrics ----
metrics.total_control_effort = trapz(sim_params.t, vecnorm(history.control_torque,2,2));
metrics = finalize_metrics(metrics, history, n_steps, sim_params.dt);

%% ---- Summary ----
display_performance_summary(metrics, history, sim_params.t);

%% ---- Plots (multi-figure + save) ----
[outdir, figdir, filelist] = generate_plots_full(history);

%% ---- Save results ----
if ~isfolder(outdir), mkdir(outdir); end
if ~isfolder(figdir), mkdir(figdir); end
stamp = datestr(now,'yyyymmdd_HHMMSS');
save(fullfile(outdir, ['simulation_' stamp '.mat']), 'history','metrics','sim_params');
fprintf('Figures exported to %s\n', figdir);

fprintf('\n==========================================\n');
fprintf('SIMULATION COMPLETE\n');
fprintf('==========================================\n');

%% ====================== LOCAL FUNCTIONS ======================
function history = store_state(history, spacecraft, index)
    history.quaternion(index,:)       = spacecraft.quaternion;
    history.angular_velocity(index,:) = spacecraft.angular_velocity;
    for i = 1:spacecraft.num_wheels
        history.wheel_speeds(index,i) = spacecraft.wheels(i).speed;
    end
end

function error_state = compute_error_state(q_ref, q_est, w_ref, w_est)
    q_err = quaternion_multiply(quaternion_conjugate(q_ref), q_est);
    if q_err(1) < 0, q_err = -q_err; end
    e_vec = q_err(2:4);
    angle_rad = 2*atan2(norm(e_vec), max(1e-12, q_err(1)));
    error_state.angle = rad2deg(angle_rad);
    error_state.vector = e_vec;
    error_state.angular_velocity = w_est - w_ref;
end

function u = apply_control_law(err, params)
    u = -params.Kp * err.vector - params.Kd * err.angular_velocity; % 1x3
end

function tau_w = allocate_to_wheels(spacecraft, body_torque_cmd)
    N = spacecraft.num_wheels;
    A = zeros(3,N);
    for i = 1:N, A(:,i) = spacecraft.wheels(i).axis.'; end
    tau_w = -pinv(A) * body_torque_cmd(:);
end

function spacecraft = propagate_dynamics(spacecraft, wheel_torques, disturbance, dt)
    q     = spacecraft.quaternion;            % 1x4
    omega = spacecraft.angular_velocity;      % 1x3
    omega_col = omega(:);
    % Quaternion kinematics
    omega_quat = [0, omega];
    q_dot = 0.5 * quaternion_multiply(omega_quat, q);
    % Wheel updates
    N = spacecraft.num_wheels; A = zeros(3,N);
    for i = 1:N
        A(:,i) = spacecraft.wheels(i).axis.';
        alpha_i = wheel_torques(i) / spacecraft.wheels(i).inertia;
        spacecraft.wheels(i).speed = spacecraft.wheels(i).speed + alpha_i * dt;
    end
    % Wheel momentum
    h_wheels = zeros(3,1);
    for i = 1:N
        h_wheels = h_wheels + spacecraft.wheels(i).axis.' * ...
                             (spacecraft.wheels(i).inertia * spacecraft.wheels(i).speed);
    end
    % Rigid body dynamics
    H_total     = spacecraft.inertia * omega_col + h_wheels;
    gyro_torque = cross(omega_col, H_total);
    applied_torque = -A * wheel_torques;
    omega_dot_col = spacecraft.inertia \ (disturbance + applied_torque - gyro_torque);
    % Integrate
    spacecraft.quaternion      = q + q_dot * dt;
    spacecraft.quaternion      = spacecraft.quaternion / norm(spacecraft.quaternion);
    spacecraft.angular_velocity = (omega_col + omega_dot_col * dt).';
end

function metrics = finalize_metrics(metrics, history, n_steps, dt)
    threshold = 2.0;
    idx = find(history.attitude_error > threshold, 1, 'last');
    if isempty(idx), metrics.settling_time = 0;
    else, metrics.settling_time = (idx-1) * dt; end
    i0 = max(1, floor(0.9*n_steps));
    metrics.steady_state_error = mean(history.attitude_error(i0:end));
    metrics.max_wheel_speed = max(max(abs(history.wheel_speeds)));
end

function display_performance_summary(metrics, history, t)
    if isempty(metrics.settling_time), metrics.settling_time = 0; end
    fprintf('========== PERFORMANCE SUMMARY ==========\n');
    fprintf('Settling time (2° threshold): %.1f seconds\n', metrics.settling_time);
    fprintf('Steady-state error: %.6f degrees\n', metrics.steady_state_error);
    fprintf('Peak attitude error: %.2f degrees\n', max(history.attitude_error));
    total_effort = trapz(t, vecnorm(history.control_torque,2,2));
    fprintf('Total control effort: %.2f N*m*s\n', total_effort);
    fprintf('Maximum wheel speed: %.1f rad/s\n', metrics.max_wheel_speed);
    fprintf('Final attitude error: %.6f degrees\n', history.attitude_error(end));
    fprintf('Final angular rates: [%.6f, %.6f, %.6f] rad/s\n', history.angular_velocity(end,:));
    fprintf('=========================================\n');
end

function [outdir, figdir, filelist] = generate_plots_full(history)
    outdir = fullfile('results');
    figdir = fullfile(outdir,'figures');
    if ~isfolder(outdir), mkdir(outdir); end
    if ~isfolder(figdir), mkdir(figdir); end
    stamp = datestr(now,'yyyymmdd_HHMMSS');
    filelist = {};

    t = history.t;
    q  = history.quaternion;
    w  = history.angular_velocity;
    ws = history.wheel_speeds;
    e  = history.attitude_error;
    u  = history.control_torque;

    % 1) Quaternion
    f1 = figure('Name','Quaternion Evolution','NumberTitle','off');
    subplot(2,2,1); plot(t, q(:,1), 'LineWidth', 1.4); grid on; xlabel('Time (s)'); ylabel('q_w'); title('Quaternion Scalar');
    subplot(2,2,2); plot(t, q(:,2), 'LineWidth', 1.4); grid on; xlabel('Time (s)'); ylabel('q_x'); title('Quaternion X');
    subplot(2,2,3); plot(t, q(:,3), 'LineWidth', 1.4); grid on; xlabel('Time (s)'); ylabel('q_y'); title('Quaternion Y');
    subplot(2,2,4); plot(t, q(:,4), 'LineWidth', 1.4); grid on; xlabel('Time (s)'); ylabel('q_z'); title('Quaternion Z');
    drawnow; exportgraphics(f1, fullfile(figdir, ['quaternion_' stamp '.png']), 'Resolution', 200);

    % 2) Angular Velocity
    f2 = figure('Name','Angular Velocity','NumberTitle','off');
    plot(t, w, 'LineWidth', 1.4); grid on; xlabel('Time (s)'); ylabel('rad/s'); title('Body Angular Rates');
    legend('\omega_x','\omega_y','\omega_z','Location','best');
    drawnow; exportgraphics(f2, fullfile(figdir, ['angular_velocity_' stamp '.png']), 'Resolution', 200);

    % 3) Attitude Error (semilog)
    f3 = figure('Name','Attitude Error','NumberTitle','off');
    semilogy(t, max(e, eps), 'LineWidth', 1.4); grid on; xlabel('Time (s)'); ylabel('deg'); title('Attitude Error vs Time');
    drawnow; exportgraphics(f3, fullfile(figdir, ['attitude_error_' stamp '.png']), 'Resolution', 200);

    % 4) Wheel Speeds
    f4 = figure('Name','Reaction Wheel Speeds','NumberTitle','off');
    plot(t, ws, 'LineWidth', 1.2); grid on; xlabel('Time (s)'); ylabel('rad/s'); title('Reaction Wheel Speeds');
    legend(arrayfun(@(i) sprintf('Wheel %d',i), 1:size(ws,2), 'UniformOutput', false),'Location','best');
    drawnow; exportgraphics(f4, fullfile(figdir, ['wheel_speeds_' stamp '.png']), 'Resolution', 200);

    % 5) Control Torque
    f5 = figure('Name','Control Torque','NumberTitle','off');
    plot(t, u, 'LineWidth', 1.2); grid on; xlabel('Time (s)'); ylabel('N*m'); title('Control Torque Components');
    legend('T_x','T_y','T_z','Location','best');
    drawnow; exportgraphics(f5, fullfile(figdir, ['control_torque_' stamp '.png']), 'Resolution', 200);
end

% ---- Quaternion utilities ----
function c = quaternion_conjugate(a)
    c = [a(1), -a(2), -a(3), -a(4)];
end
function c = quaternion_multiply(a, b)
    w1=a(1); x1=a(2); y1=a(3); z1=a(4);
    w2=b(1); x2=b(2); y2=b(3); z2=b(4);
    c = [ w1*w2 - x1*x2 - y1*y2 - z1*z2, ...
          w1*x2 + x1*w2 + y1*z2 - z1*y2, ...
          w1*y2 - x1*z2 + y1*w2 + z1*x2, ...
          w1*z2 + x1*y2 - y1*x2 + z1*w2 ];
end
