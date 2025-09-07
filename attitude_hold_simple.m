%% Simplified Spacecraft Attitude Control Simulation
% This is a simplified, working version to ensure basic functionality
clear all; close all; clc;

fprintf('==========================================\n');
fprintf('SPACECRAFT ATTITUDE CONTROL SIMULATION\n');
fprintf('(Simplified Version)\n');
fprintf('==========================================\n\n');

%% Spacecraft Parameters
% Inertia matrix (kg*m^2)
I_spacecraft = diag([6400, 4730, 8160]);

% Initial conditions
q_init = [sqrt(2)/2, 0, 0, sqrt(2)/2];  % Initial quaternion
omega_init = [0.01, 0.01, 0.01];        % Initial angular velocity (rad/s)

% Reaction wheel parameters
num_wheels = 4;
wheel_inertia = 2.0;  % kg*m^2
wheel_axes = [1, 0, 0;      % Wheel 1: X-axis
              0, 1, 0;      % Wheel 2: Y-axis
              0, 0, 1;      % Wheel 3: Z-axis
              1, 1, 1]/sqrt(3);  % Wheel 4: Diagonal

% Controller gains
Kp = 10 * eye(3);
Kd = 450 * eye(3);

%% Simulation Parameters
t_final = 1000;  % seconds
dt = 0.1;        % seconds
t = 0:dt:t_final;
n_steps = length(t);

%% Initialize State Variables
% Current state
q = q_init;
omega = omega_init;
wheel_speeds = zeros(1, num_wheels);

% Target state
q_target = [1, 0, 0, 0];  % Identity quaternion
omega_target = [0, 0, 0];

% History arrays
q_history = zeros(n_steps, 4);
omega_history = zeros(n_steps, 3);
wheel_speed_history = zeros(n_steps, num_wheels);
error_history = zeros(n_steps, 1);
torque_history = zeros(n_steps, 3);

%% Main Simulation Loop
fprintf('Running simulation...\n');
fprintf('Time: ');

for k = 1:n_steps
    % Store current state
    q_history(k,:) = q;
    omega_history(k,:) = omega;
    wheel_speed_history(k,:) = wheel_speeds;
    
    % === CONTROL ===
    % Compute quaternion error
    q_error = quat_multiply(q, quat_conjugate(q_target));
    if q_error(1) < 0
        q_error = -q_error;
    end
    
    % Extract error vector
    error_vector = q_error(2:4);
    error_history(k) = 2 * acos(min(1, abs(q_error(1)))) * 180/pi;
    
    % PD control law
    omega_error = omega - omega_target;
    control_torque = -Kp * error_vector(:) - Kd * omega_error(:);
    torque_history(k,:) = control_torque';
    
    % === CONTROL ALLOCATION ===
    % Build wheel matrix (3 x num_wheels)
    A = wheel_axes';
    
    % Current wheel momentum
    h_wheels = wheel_inertia * wheel_speeds(:);
    
    % Gyroscopic compensation
    omega_cross = [0, -omega(3), omega(2);
                   omega(3), 0, -omega(1);
                   -omega(2), omega(1), 0];
    h_total_wheels = A * h_wheels;
    gyro_comp = A' * pinv(A * A') * omega_cross * h_total_wheels;
    
    % Wheel torque commands
    wheel_torques = -pinv(A) * control_torque - gyro_comp;
    
    % === DYNAMICS ===
    if k < n_steps
        % Quaternion kinematics
        omega_quat = [0, omega];
        q_dot = 0.5 * quat_multiply(omega_quat, q);
        
        % Update wheel speeds
        wheel_accelerations = wheel_torques / wheel_inertia;
        wheel_speeds = wheel_speeds + wheel_accelerations' * dt;
        
        % Angular momentum
        h_wheels_new = wheel_inertia * wheel_speeds(:);
        h_total_wheels_new = A * h_wheels_new;
        
        % Euler's equation
        H_total = I_spacecraft * omega(:) + h_total_wheels_new;
        applied_torque = -A * wheel_torques;
        gyro_torque = cross(omega(:), H_total);
        
        omega_dot = I_spacecraft \ (applied_torque - gyro_torque);
        
        % Integrate
        q = q + q_dot * dt;
        q = q / norm(q);  % Normalize
        omega = omega + omega_dot' * dt;
    end
    
    % Progress display
    if mod(k, floor(n_steps/10)) == 0
        fprintf('%.0f%% ', 100*k/n_steps);
    end
end

fprintf('\nSimulation complete!\n\n');

%% Display Results
fprintf('========== RESULTS ==========\n');
fprintf('Initial attitude error: %.3f degrees\n', error_history(1));
fprintf('Final attitude error: %.6f degrees\n', error_history(end));
fprintf('Final angular rates: [%.6f, %.6f, %.6f] rad/s\n', omega_history(end,:));
fprintf('Final wheel speeds: [%.3f, %.3f, %.3f, %.3f] rad/s\n', wheel_speed_history(end,:));
fprintf('=============================\n\n');

%% Plot Results
% Figure 1: Quaternion
figure('Name', 'Quaternion Evolution');
subplot(2,2,1);
plot(t, q_history(:,1), 'LineWidth', 2);
grid on; xlabel('Time [s]'); ylabel('q_w');
title('Quaternion Scalar');

subplot(2,2,2);
plot(t, q_history(:,2), 'LineWidth', 2);
grid on; xlabel('Time [s]'); ylabel('q_x');
title('Quaternion X');

subplot(2,2,3);
plot(t, q_history(:,3), 'LineWidth', 2);
grid on; xlabel('Time [s]'); ylabel('q_y');
title('Quaternion Y');

subplot(2,2,4);
plot(t, q_history(:,4), 'LineWidth', 2);
grid on; xlabel('Time [s]'); ylabel('q_z');
title('Quaternion Z');

% Figure 2: Angular Velocity
figure('Name', 'Angular Velocity');
plot(t, omega_history(:,1), 'b-', 'LineWidth', 2); hold on;
plot(t, omega_history(:,2), 'r-', 'LineWidth', 2);
plot(t, omega_history(:,3), 'g-', 'LineWidth', 2);
grid on; xlabel('Time [s]'); ylabel('Angular Velocity [rad/s]');
legend('\omega_x', '\omega_y', '\omega_z');
title('Angular Velocity Components');

% Figure 3: Attitude Error
figure('Name', 'Attitude Error');
semilogy(t, error_history, 'LineWidth', 2);
grid on; xlabel('Time [s]'); ylabel('Error [degrees]');
title('Attitude Error vs Time');

% Figure 4: Wheel Speeds
figure('Name', 'Reaction Wheel Speeds');
plot(t, wheel_speed_history, 'LineWidth', 2);
grid on; xlabel('Time [s]'); ylabel('Wheel Speed [rad/s]');
legend('Wheel 1', 'Wheel 2', 'Wheel 3', 'Wheel 4');
title('Reaction Wheel Speeds');

% Figure 5: Control Torque
figure('Name', 'Control Torque');
plot(t, torque_history, 'LineWidth', 2);
grid on; xlabel('Time [s]'); ylabel('Torque [N*m]');
legend('T_x', 'T_y', 'T_z');
title('Control Torque Components');

%% Helper Functions
function q_out = quat_multiply(q1, q2)
    % Quaternion multiplication
    w1 = q1(1); x1 = q1(2); y1 = q1(3); z1 = q1(4);
    w2 = q2(1); x2 = q2(2); y2 = q2(3); z2 = q2(4);
    
    q_out = [w1*w2 - x1*x2 - y1*y2 - z1*z2, ...
             w1*x2 + x1*w2 + y1*z2 - z1*y2, ...
             w1*y2 - x1*z2 + y1*w2 + z1*x2, ...
             w1*z2 + x1*y2 - y1*x2 + z1*w2];
end

function q_conj = quat_conjugate(q)
    % Quaternion conjugate
    q_conj = [q(1), -q(2), -q(3), -q(4)];
end