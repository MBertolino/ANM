clear all; close all;

% Choose order of SBP
ordning = 4;

% Setup time
t_start = 0;
t_end = 1.8;
t = t_start;

% Setup space
x_l = -1;
x_r = 1;
x0 = 0;
L_x = x_r - x_l;
y_b = -1;
y_t = 1;
y0 = 0;
L_y = y_t - y_b;

% Setup grid
m = 51;
x = linspace(x_l, x_r, m); % Discrete x-values
y = linspace(y_b, y_t, m);
h = L_x/(2*m-1);
dt = 0.1*h;
N_iter = floor(t_end/dt);

% Setup system
A = [0 0 0; 0 0 -1; 0 -1 0];
B = [0 1 0; 1 0 0; 0 0 0];
tau_w = [0; -1; -2]; % Penalty parameter
tau_e = [0; -1; 2];
tau_n = [-2; -1; 0];
tau_s = [2; -1; 0];
e_1u = [1 0 0]; % Choose variable
e_2u = [0 1 0];
e_3u = [0 0 1];

% Setup Gaussian width
rr = 0.1;

% Load operators
Val_operator_ANM;
eye_m = eye(m);
D_x = kron(D1, eye_m);
D_y = kron(eye_m, D1);
H_x = kron(H, eye_m);
H_y = kron(eye_m, H);
HI_x = inv(H_x);
HI_y = inv(H_y);

% Variables
e_w = kron(e_1, eye_m);
e_e = kron(e_m, eye_m);
e_s = kron(eye_m, e_1);
e_n = kron(eye_m, e_m);

% SBP = -SAT approximation
SAT_w = kron(tau_w, HI_x)*e_w*(kron(e_2u, e_w'));
SAT_e = kron(tau_e, HI_x)*e_e*(kron(e_2u, e_e'));
SAT_s = kron(tau_s, HI_y)*e_s*(kron(e_2u, e_s'));
SAT_n = kron(tau_n, HI_y)*e_n*(kron(e_2u, e_n'));
PP = kron(A, D_x) + kron(B, D_y) + SAT_w + SAT_e + SAT_s + SAT_n;
P = dt*sparse(PP);

% Initial position
V = [zeros(m*m, 1); theta_init(x, x0, y, y0, rr); zeros(m*m, 1)];

% Pre-allocate RK-vectors and errors
temp = zeros(3*m*m, 1);    % Temporary vector in RK4
w1 = zeros(3*m*m, 1);      % Step 1 vector in RK4
w2 = zeros(3*m*m, 1);      % Step 2 vector in RK4
w3 = zeros(3*m*m, 1);      % Step 3 vector in RK4
w4 = zeros(3*m*m, 1);      % Step 4 vector in RK4

% RK4
for k = 1:N_iter
    w1 = P*V;
    temp = V + w1*0.5;
    
    w2 = P*temp;
    temp = V + w2*0.5;
    
    w3 = P*temp;
    temp = V + w3;
    
    w4 = P*temp;
    
    V = V + (w1 + 2*w2 + 2*w3 + w4)/6;
    Vmat_1 = vec2mat(V(1:m*m), m);
    Vmat_2 = vec2mat(V(m*m+1:2*m*m), m);
    Vmat_3 = vec2mat(V(2*m*m+1:3*m*m), m);
    
    t = t + dt;
    
    % Plot
    surf(x, y, Vmat_1)
    xlim([-1 1])
    ylim([-1 1])
    zlim([-1 1])
    pause(0.001)
end

% Help functions
function theta = theta_init(x, x0, y, y0, rr)
m = length(x);
n = length(y);
theta = zeros(m*n, 1);
for (ix = 0:(m-1))
    for (iy = 1:n)
        theta(ix*m + iy, 1) = exp(-((x(ix+1) - x0)/rr).^2 - ((y(iy) - y0)/rr).^2)';
    end
end
end
