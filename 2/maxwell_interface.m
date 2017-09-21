clear all; close all;

% Choose order and set grid refine_muents
% grid_ref = [1 2 3 4];
grid_ref = 8;

% Repeat for different orders of SBP
% order = [2 4 6 10];
ordning = 4;

% Setup time
t_start = 0;
t_end = 1.8;

% Setup space
x_l = -1;
x_i = 0;
x_r = 1;
x_0 = x_l*0.5;
L = x_r - x_l;

% Setup syste_mu
A = [0 1; 1 0];
[S, Lambda] = eig(A);
Lambda_pos = (Lambda + abs(Lambda))*0.5;
Lambda_neg = (Lambda - abs(Lambda))*0.5;
A_pos = S*Lambda_pos/S;
A_neg = S*Lambda_neg/S;

epsilon_l = 2;
epsilon_r = 1;
eta_l = sqrt(epsilon_l); % Refraction
eta_r = sqrt(epsilon_r);
T_exact = 2*eta_l./(eta_l + eta_r); % Transmission
R_exact = (eta_l - eta_r)./(eta_l + eta_r); % Reflection
C_l = [epsilon_l 0; 0 1];
C_r = [epsilon_r 0; 0 1];
C_lI = inv(C_l);
C_rI = inv(C_r);

sigma_l = 1;
sigma_r = -1;
tau_l = [-1; 1]; % Penalty parameter
tau_r = [-1; -1];
e_1u = [1 0]; % Choose variable
e_ku = [0 1];

% Setup exact solution
rr = 0.1;          % Width of Gaussian

% Repeat for two grid refinements to check convergence
for j = 1:length(grid_ref) % Grid refinements for convergence
    
    % Start from 0
    t = t_start;
    
    % Setup grid
    m = grid_ref(j)*31;
    x = linspace(x_l, x_r, m); % Discrete x-values
    x_L = linspace(x_l, 0, m);
    x_R = linspace(0, x_r, m);
    h = L/(2*m-1);
    dt = 0.1*h;
%     N_iter = floor(t_end/dt);
    N_iter = round(((0.5*eta_l)/dt) + (0.5*eta_r)/dt); % N_iter s.t. waves stop at 0.5
    
    % Load operators
    Val_operator_ANM;
    
    % SBP = -SAT approximation
    PP_l = kron(C_l\A, D1) + sigma_r*kron(C_l\A_pos, HI*e_m*e_m') + ...
        + kron(C_l\tau_l*e_1u, HI*e_1*e_1');
    PP_r = kron(C_r\A, D1) + sigma_l*kron(C_r\A_neg, HI*e_1*e_1') + ...
        + kron(C_r\tau_r*e_1u, HI*e_m*e_m');
    
    P_l = dt*sparse(PP_l);
    P_r = dt*sparse(PP_r);
    
    G_l = dt*sigma_r*sparse(-kron(C_l\A_pos, HI))*kron(eye(2), e_m*e_1');    % Penalty data
    G_r = dt*sigma_l*sparse(-kron(C_r\A_neg, HI))*kron(eye(2), e_1*e_m');
    
    % Initial position
    V_l = [theta_2(x_L, x_0, 0, rr) - theta_1(x_L, x_0, 0, rr); theta_2(x_L, x_0, 0, rr) + theta_1(x_L, x_0, 0, rr)];
    V_r = zeros(2*m, 1);
    
    % Exact solution
    %         exact = [theta_1(x_L, x_0, L-t, rr) - theta_2(x_L, x_0, L-t, rr); theta_1(x_L, x_0, L-t, rr) + theta_2(x_L, x_0, L-t, rr)];
    
    % Pre-allocate RK-vectors and errors
    temp_l = zeros(2*m, 1);    % Te_muporary vector in RK4
    w1_l = zeros(2*m, 1);      % Step 1 vector in RK4
    w2_l = zeros(2*m, 1);      % Step 2 vector in RK4
    w3_l = zeros(2*m, 1);      % Step 3 vector in RK4
    w4_l = zeros(2*m, 1);      % Step 4 vector in RK4
    temp_r = zeros(2*m, 1);    % Te_muporary vector in RK4
    w1_r = zeros(2*m, 1);      % Step 1 vector in RK4
    w2_r = zeros(2*m, 1);      % Step 2 vector in RK4
    w3_r = zeros(2*m, 1);      % Step 3 vector in RK4
    w4_r = zeros(2*m, 1);      % Step 4 vector in RK4
    
    % Setup video
    vidObj = VideoWriter('Test_1D.avi');
    open(vidObj);
    update_movie = 10; % How often to update movie
    
    % RK4 on left and right side
    for k = 1:N_iter
        w1_l = P_l*V_l + G_l*V_r;
        w1_r = P_r*V_r + G_r*V_l;
        temp_l = V_l + w1_l*0.5;
        temp_r = V_r + w1_r*0.5;
        
        w2_l = P_l*temp_l + G_l*temp_r;
        w2_r = P_r*temp_r + G_r*temp_l;
        temp_l = V_l + w2_l*0.5;
        temp_r = V_r + w2_r*0.5;
        
        w3_l = P_l*temp_l + G_l*temp_r;
        w3_r = P_r*temp_r + G_r*temp_l;
        temp_l = V_l + w3_l;
        temp_r = V_r + w3_r;
        
        w4_l = P_l*temp_l + G_l*temp_r;
        w4_r = P_r*temp_r + G_r*temp_l;
        
        V_l = V_l + (w1_l + 2*w2_l + 2*w3_l + w4_l)/6;
        V_r = V_r + (w1_r + 2*w2_r + 2*w3_r + w4_r)/6;
        
        t = t + dt;

        % Plot
%         if (mod(k, update_movie) == 0)
%             plot(x_L, V_l(1:m), 'b', x_L, V_l(m+1:end), 'r', ...
%                 x_R, V_r(1:m), 'b', x_R, V_r(m+1:end), 'r')
%             ylim([-2 2])
%             currFrame = getframe;
%             writeVideo(vidObj, currFrame);
%         end
    end
    close(vidObj)
end

% Calculate transmissions
E_I = 1; % Incident wave
E_tx = abs(min(V_r(1:m))); % Transmitted wave
E_rx = abs(min(V_l(1:m)));
T = E_tx/E_I
R = E_rx/E_I

err_T = abs(T - T_exact);
err_R = abs(R - R_exact);
err = 0.5*(err_T + err_R)

% Help functions
function theta = theta_1(x, x0, t, rr)
theta = exp(-((x - x0 - t)/rr).^2)';
end

function theta = theta_2(x, x0, t, rr)
theta = -exp(-((x - x0 + t)/rr).^2)';
end
