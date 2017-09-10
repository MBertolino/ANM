clear all; close all;

% Setup space
m = 31;
x_l = -1;
x_r = 1;
L = x_r - x_l;
x = linspace(x_l, x_r, m); % Discrete x-values
h = L/(m-1);

% Setup time
t_start = 0;
t = t_start;
t_end = 1.8;                
dt = 0.1*h;
N_iter = t_end/dt;

% Setup system
A = [0 1; 1 0];
[S, Lambda] = eig(A);
Lambda_pos = (Lambda + abs(Lambda))*0.5;
Lambda_neg = (Lambda - abs(Lambda))*0.5;
A_pos = S*Lambda_pos/S;
A_neg = S*Lambda_neg/S;
% tau_l = [-1; 1]; % Penalty parameter
% tau_r = [-1; -1];
e_1u = [1 0]; % Choose variable
e_ku = [0 1];

% Setup operators
ordning = 6;
Val_operator_ANM;

% SBP = -SAT approximation for Dirichlet
% PP = kron(A, D1) + kron(tau_l, HI)*e_1*kron(e_1u, e_1') + ...
%     kron(tau_r, HI)*e_m*kron(e_1u, e_m');
% P = dt*sparse(PP);
% G_l = dt*sparse(-kron(tau_l, HI*e_1));    % Penalty data
% G_r = dt*sparse(-kron(tau_r, HI*e_m));

% SBP = -SAT approximation for characteristic
PP = kron(A, D1) + kron(A_neg, HI*e_1*e_1') - kron(A_pos, HI*e_m*e_m');
P = dt*sparse(PP);

% Setup exact solution
rr = 0.1;          % Width of Gaussian
gl = 0;
gr = 0;
V = [theta_2(x, 0, rr) - theta_1(x, 0, rr); theta_2(x, 0, rr) + theta_1(x, 0, rr)];
exact = [theta_1(x, L-t, rr) - theta_2(x, L-t, rr); theta_1(x, L-t, rr) + theta_2(x, L-t, rr)];

% Pre-allocate RK-vectors and errors
temp = zeros(2*m, 1);    % Temporary vector in RK4
w1 = zeros(2*m, 1);      % Step 1 vector in RK4
w2 = zeros(2*m, 1);      % Step 2 vector in RK4
w3 = zeros(2*m, 1);      % Step 3 vector in RK4
w4 = zeros(2*m, 1);      % Step 4 vector in RK4
err_1 = zeros(m, 1);
err_2 = zeros(m, 1);

% Setup video
vidObj = VideoWriter('Test_1D.avi');
open(vidObj);
update_movie = 10; % How often to update movie

for k = 1:N_iter
    %%%--------Start R-K without dissipation-------------------------------
    w1 = P*V;% + G_l*gl + G_r*gr;
    temp = V + w1/2;
    
    %     gl = exp(-((x_l - (t+dt/2) - x_0)/rr).^2)';
    %     gr = exp(-((x_r - (t+dt/2) - x_0)/rr).^2)';
    w2 = P*temp; % + G_l*gl + G_r*gr;
    temp = V + w2/2;
    
    w3 = P*temp; % + G_l*gl + G_r*gr;
    temp = V + w3;
    
    %     gl = exp(-((x_l - (t+dt) - x_0)/rr).^2)';
    %     gr = exp(-((x_l - (t+dt) - x_0)/rr).^2)';
    w4 = P*temp;% + G_l*gl + G_r*gr;
    
    V = V + (w1 + 2*w2 + 2*w3 + w4)/6;
    t = t + dt;
    exact = [theta_1(x, L-t, rr) - theta_2(x, L-t, rr); ...
        theta_1(x, L-t, rr) + theta_2(x, L-t, rr)];
    
    % Plot
%     if mod(k, update_movie) == 0
%         plot(x, V(1:m), x, V(m+1:end), 'r')
%         ylim([-2 2])
%         currFrame = getframe;
%         writeVideo(vidObj, currFrame);
%     end
    
    % Calculate error
    err_1(k) = sqrt(h)*norm(V(1:m) - exact(1:m));
    err_2(k) = sqrt(h)*norm(V(m+1:end) - exact(m+1:end));
end
close(vidObj)

% Plot eigenvalues
figure()
plot(real(eig(PP)), imag(eig(PP)), '*')
xlabel('Real axis')
ylabel('Imaginary axis')

% Help functions
function theta = theta_1(x, t, rr)
    theta = exp(-((x - t)/rr).^2)';
end

function theta = theta_2(x, t, rr)
    theta = -exp(-((x + t)/rr).^2)';
end
