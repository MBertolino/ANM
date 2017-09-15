
clear all;
close all;

% Test different orders of SBP
ordning_all = [2 4 6 10];
ordning = 2;
m = 40;

rr = 0.1;          % Width of Gaussian
x_0 = -0.5;           % Initial position of Gaussian
x_l = -1;
x_r = 1;    % The boundaries of the domain
L = x_r - x_l;
e2 = [0 1];

xL = linspace(x_l, 0, m);
xR = linspace(0, x_r, m);
t = 0;
A = [0 1; 1 0];
[S, E] = eig(A);
A_p = 0.5*S'*(E+abs(E))*S;
A_m = 0.5*S'*(E-abs(E))*S;
CR = [0.5 0;0 1];
CRI = inv(CR);
CL = [2 0; 0 1];
CLI = inv(CL);
VL = zeros(1, 2*m);
VR = zeros(1, 2*m);
e1 = [1 0];

% Initial condition
%VL(1:m) = exp(-((xL-t)/rr).^2);
%VL(m+1:end) = -exp(-((xL+t)/rr).^2)*0;
%VR(1:m) = exp(-((xR-t)/rr).^2);
%VR(m+1:end) = -exp(-((xR+t)/rr).^2);

%1.0
t_1 = 1.0;           % End time

tau_l = [0; 1];          % Penalty parameter
tau_r = [0; -1];
sigma_l = 1;
sigma_r = -1;

h = L/(2*m - 1);
Val_operator_ANM;
t = 0;
ff = 0.001;
dt = ff*h;                 % Time-step used
max_itter = floor(t_1/dt);  % Number of steps in time
gl = 0;
gr = 0;
PPL = kron(CLI*A, D1) + kron(CLI*tau_l, HI*e_1)*kron(e1, e_1') + ...
    sigma_r*kron(CLI*A_p, HI)*kron(eye(2), e_m*e_m'); % SBP=-SAT approximation
GL_R = -dt*sigma_r*kron(CLI*A_p, HI)*kron(eye(2), e_m*e_1');

PPR = kron(CRI*A, D1) + kron(CRI*tau_r, HI*e_m)*kron(e1, e_m') + ...
    sigma_l*kron(CRI*A_m, HI)*kron(eye(2), e_1*e_1');
GR_L = -dt*sigma_l*kron(CRI*A_m, HI)*kron(eye(2), e_1*e_m');

PL = dt*sparse(PPL);
PR = dt*sparse(PPR);

tempL = zeros(m, 1);    % Temporary vector in RK4
tempR = zeros(m, 1);

wL1 = zeros(m, 1);      % Step 1 vector in RK4
wL2 = zeros(m, 1);      % Step 2 vector in RK4
wL3 = zeros(m, 1);      % Step 3 vector in RK4
wL4 = zeros(m, 1);      % Step 4 vector in RK4
wR1 = zeros(m, 1);      % Step 1 vector in RK4
wR2 = zeros(m, 1);      % Step 2 vector in RK4
wR3 = zeros(m, 1);      % Step 3 vector in RK4
wR4 = zeros(m, 1);      % Step 4 vector in RK4

felet=zeros(max_itter+1,1); % Error vector in time

exactfun=@(x,t)[exp(-((x-L+t-x_0)/rr).^2)+exp(-((x+L-t-x_0)/rr).^2), exp(-((x-L+t-x_0)/rr).^2)-exp(-((x+L-t-x_0)/rr).^2)]';    %Exact solution
exact=exactfun(xL,t);
VL=[-exp(-((xL-x_0)/rr).^2)-exp(-((xL-x_0)/rr).^2), exp(-((xL-x_0)/rr).^2)-exp(-((xL-x_0)/rr).^2)]';% Numerical solution
VR=[-exp(-((xR-x_0)/rr).^2)-exp(-((xR-x_0)/rr).^2), exp(-((xR-x_0)/rr).^2)-exp(-((xR-x_0)/rr).^2)]';% Numerical solution


felet(1)=sqrt(h)*norm(VL-exact);
%gl=exp(-((x_l-t-x_0)/rr).^2)';
xn=xL;
%x=[linspace(x_l,x_r,m),linspace(x_l,x_r,m)];


% Setup video
vidObj = VideoWriter('Test_1D.avi');
open(vidObj);
update_movie = 1/ff; % How often to update movie

for nr_itter=1:max_itter
    
    %%%--------Start R-K utan dissipation----------------------------------------
    wL1 = PL*VL + GL_R*VR;
    wR1 = PR*VR + GR_L*VL;
    tempL = VL + wL1*0.5;
    tempR = VR + wR1*0.5;
    
    wL2 = PL*tempL + GL_R*tempR;
    wR2 = PR*tempR + GR_L*tempL; 
    tempL = VL + wL2*0.5;
    tempR = VR + wR2*0.5;
    
    wL3 = PL*tempL + GL_R*tempR;
    wR3 = PR*tempR + GR_L*tempL;
    tempL = VL + wL3;
    
    % kom ihåg att uppdatera gl och gr
    wL4 = PL*tempL + GL_R*tempR;
    wR4 = PR*tempR + GR_L*tempL;
    
    VL = VL + (wL1 + 2*wL2 + 2*wL3 + wL4)/6;
    VR = VR + (wR1 + 2*wR2 + 2*wR3 + wR4)/6;
    
    t = t + dt;
    if(abs(t-0.35)<=0.9*dt)
        mL1 = round(0.72/h);
        mL2 = round(0.76/h);
        Emax = abs(min(VL(mL1:mL2)));
    end
    
    exact=exactfun(xL,t);
    felet(nr_itter+1)=sqrt(h)*norm(VL-exact);
    if mod(nr_itter, update_movie) == 0
        plot(xL, VL(1:m),'b', xL, VL(m+1:end), 'r',xR, VR(1:m),'b', xR, VR(m+1:end), 'r')
        ylim([-2 2])
        currFrame = getframe;
        writeVideo(vidObj, currFrame);
    end
    
end
close(vidObj)

plot(real(eig(PPR)), imag(eig(PPR)),'*');
xlabel('Real');
ylabel('Imaginary');

etaL = sqrt(CL(1,1));
etaR = sqrt(CR(1,1));
Texact = 2*etaL/(etaL+etaR);
Rexact = (etaL-etaR)/(etaL+etaR);
mL1 = ceil(0.76/h);
mL2 = ceil(0.82/h);
mR1 = ceil(0.36/h);
mR2 = ceil(0.45/h);
ELmax = abs(min(VL(mL1:mL2)));
ERmax=abs(min(VR(mR1:mR2)));

T = ERmax/Emax;
R = ELmax/Emax;
fel = abs(T-Texact);
fel2 = abs(R-Rexact);
(fel+fel2)/2


