clear all; close all;

% Choose order and set grid refine_muents
% grid_ref = [1 2 3 4];
grid_ref = 4;

% Repeat for different orders of SBP
% order = [2 4 6 10];
order = 4;
for iOrder = 1:length(order)
    ordning = order(iOrder);
    
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
    eta_l = 2; % Refraction
    eta_r = 1;
    T = 2*eta_l./(eta_l + eta_r); % Transmission
    R = (eta_l - eta_r)./(eta_l + eta_r); % Reflection
    A = [0 1; 1 0];
    [S, Lambda] = eig(A);
    Lambda_pos = (Lambda + abs(Lambda))*0.5;
    Lambda_neg = (Lambda - abs(Lambda))*0.5;
    A_pos = S*Lambda_pos/S;
    A_neg = S*Lambda_neg/S;
    C_l = [2 0; 0 1];
    C_lI = inv(C_l);
    C_r = [1 0; 0 1];
    C_rI = inv(C_r);
    sigma_l = 1;
    sigma_r = -1;
    tau_l = [-1; 1]; % Penalty parameter
    tau_r = [-1; -1];
    e_1u = [1 0]; % Choose variable
    e_ku = [0 1];
    
    % Setup exact solution
    rr = 0.1;          % Width of Gaussian
    
    % Repeat for two grid refine_muents to check convergence
    for j = 1:length(grid_ref) % Grid refine_muents for convergence
        
        % Start from 0
        t = t_start;
        
        % Setup grid
        m = grid_ref(j)*31;
        x = linspace(x_l, x_r, m); % Discrete x-values
        x_L = linspace(x_l, 0, m);
        x_R = linspace(0, x_r, m);
        h = L/(2*m-1);
        dt = 0.1*h;
        N_iter = t_end/dt;
        %     err = zeros(int16(N_iter-1), length(grid_ref));
        
        % Load operators
        Val_operator_ANM;
        
        % SBP = -SAT approximation
        PP_l = kron(C_l\A, D1) + sigma_r*kron(C_l\A_pos, HI*e_m*e_m') + ...
            + kron(C_l\tau_l*e_1u, HI*e_1*e_1');
        PP_r = kron(C_r\A, D1) + sigma_l*kron(C_r\A_neg, HI*e_1*e_1') + ...
            + kron(C_r\tau_r*e_1u, HI*e_m*e_m');
        
        P_l = dt*sparse(PP_l);
        P_r = dt*sparse(PP_r);
        
        G_l = dt*sigma_r*sparse(-kron(C_lI*A_pos, HI))*kron(eye(2), e_m*e_1');    % Penalty data
        G_r = dt*sigma_l*sparse(-kron(C_rI*A_neg, HI))*kron(eye(2), e_1*e_m');
        
        % Analytical solution
        V_l = [theta_2(x_L, x_0, 0, rr) - theta_1(x_L, x_0, 0, rr); theta_2(x_L, x_0, 0, rr) + theta_1(x_L, x_0, 0, rr)];
        V_r = [theta_2(x_R, x_0, 0, rr) - theta_1(x_R, x_0, 0, rr); theta_2(x_R, x_0, 0, rr) + theta_1(x_R, x_0, 0, rr)];
        
        % Pre-allocate RK-vectors and errors
        temp_l = zeros(m, 1);    % Te_muporary vector in RK4
        w1_l = zeros(m, 1);      % Step 1 vector in RK4
        w2_l = zeros(m, 1);      % Step 2 vector in RK4
        w3_l = zeros(m, 1);      % Step 3 vector in RK4
        w4_l = zeros(m, 1);      % Step 4 vector in RK4
        temp_r = zeros(m, 1);    % Te_muporary vector in RK4
        w1_r = zeros(m, 1);      % Step 1 vector in RK4
        w2_r = zeros(m, 1);      % Step 2 vector in RK4
        w3_r = zeros(m, 1);      % Step 3 vector in RK4
        w4_r = zeros(m, 1);      % Step 4 vector in RK4
        
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
            %             exact = [theta_1(x, 0, L-t, rr) - theta_2(x, 0, L-t, rr); ...
            %                 theta_1(x, 0, L-t, rr) + theta_2(x, 0, L-t, rr)];
            
            % Plot
            if mod(k, update_movie) == 0
                plot(x_L, V_l(1:m), 'b', x_L, V_l(m+1:end), 'r', ...
                    x_R, V_r(1:m), 'b', x_R, V_r(m+1:end), 'r') %, x, exact(1:m), x, exact(m+1:end))
                ylim([-2 2])
                currFrame = getframe;
                writeVideo(vidObj, currFrame);
            end
            
            % Calculate error
            %         err(k, j) = sqrt(h)*norm(V - exact);
        end
        close(vidObj)
        %         err(j) = sqrt(h)*norm(V_l - exact);
        
        %         if j > 1
        %             q(iOrder, j-1) = log10(err(j)/err(j-1))/log10(grid_ref(j)/grid_ref(j-1));
        %         end
    end
    
    % Check convergence
    % conv = log10(err(:, 2)./err(:, 1))/log10(grid_ref(2)/grid_ref(1));
end

% % Plot convergence
% figure()
% plot(1:N_iter, conv)
% xlabel('Iterations')
% ylabel('q')

% Plot eigenvalues
figure()
plot(real(eig(PP_l)), imag(eig(PP_l)), '*')
xlabel('Real axis')
ylabel('Imaginary axis')

% Help functions
function theta = theta_1(x, x0, t, rr)
theta = exp(-((x - x0 - t)/rr).^2)';
end

function theta = theta_2(x, x0, t, rr)
theta = -exp(-((x - x0 + t)/rr).^2)';
end
