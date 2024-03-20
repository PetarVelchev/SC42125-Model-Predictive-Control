%% SC42125 Model Predictive Control - Project 
% Authors: Floris van Leeuwen, Petar Velchev 
% Department for Systems & Control 
% Delft Univirsity of Technology, The Netherlands 

%% Clean and load 
clc;
clear all;
close all;
%load vehicle parameters from m file
vehicle_params;

%% Initial conditions 
Vx  =  75/3.6;              % vehicle speed, km/h
SWA = 30;                   % steering wheel angle, deg

%% Prepare cornering stiffness tire model
Ct  = Ca_f + Ca_r;
Cs  = l_f * Ca_f - l_r * Ca_r;
Cq2 = l_f^2 * Ca_f + l_r^2 * Ca_r;

%% Continuous system dynamics 

%Simple model:
% % state vector x = [v r]
% % control vector u = [delta]
% % output vector y = [ay r beta v]
% % delta - road steering angle, rad
% % ay - lateral acceleration, m/s^2
% % r - yaw rate, rad/s
% % v - lateral velocity, m/s
% % beta - sideslip angle, rad
% 
% A_ss = -[Ct / (m * Vx), Cs / (m * Vx) + Vx; Cs / (Izz * Vx), Cq2 / (Izz * Vx)];
% B_ss = [Ca_f / m; l_f * Ca_f / Izz];
% C_ss = [-Ct / (m * Vx), -Cs / (m * Vx); 0, 1; 1/Vx, 0; 1 0];
% D_ss = [Ca_f / m; 0; 0; 0];

%More complex model:
% state vector x = [y vy yaw r]
% control vector u = [delta]
% output vector y = [y ay r vy]
% delta - road steering angle, rad
% ay - lateral acceleration, m/s^2
% r - yaw rate, rad/s
% vy - lateral velocity, m/s
% vx - longitudinal velocity, m/s

A_ss = -[0, 1, 0, 0; 
    0, -Ct / (m * Vx), 0,-Vx - Cs/(m * Vx);
    0, 0, 0, 1;
    0, -Cs/(Izz * Vx), 0, -Cq2/(Izz * Vx)];
B_ss = [0;
    Ca_f / m; 
    0;
    l_f * Ca_f / Izz];
C_ss = [1, 0, 0, 0;
    0, -Ct / (m * Vx), 0,-Cs / (m * Vx); 
    0, 0, 0, 1; 
    0, 1, 0, 0];
D_ss = [0; Ca_f / m; 0; 0];

% System dimensions


% Define transfer function
sys_cont = ss(A_ss,B_ss,C_ss,D_ss);
systf_cont = tf(sys_cont);

%% Discrete system dynamics 
T =5; %Simulation Time
Ts = 0.1; %Sampling Time
t_span = 0:Ts:T; %span vector
%x0 = [0, 0, 0, 0]';
simulation_steps = T/Ts; 

sys_disc = c2d(sys_cont, Ts);

% LTI setup
LTI.A = sys_disc.A;
LTI.B = sys_disc.B;
LTI.C = sys_disc.C;
LTI.D = sys_disc.D;

LTI.yref=[2; 0; 0; 0];             
LTI.x0=[0; 0; 0; 0];

%Definition of system dimension
dim.nx=4;     %state dimension
dim.nu=1;     %input dimension
dim.ny=4;     %output dimension
dim.nd=0;     %disturbance dimension
dim.N=10;      %horizon

%Definition of quadratic cost function
weight.Q=[10 0 0 0; 0 1 0 0; 0 0 3 0; 0 0 0 1];                %weight on output
weight.R=eye(dim.nu);                          %weight on input
weight.P=dare(LTI.A,LTI.B,weight.Q,weight.R);  %terminal cost

predmod=predmodgen(LTI,dim);  
[H,h]=costgen(predmod,weight,dim); 
%L=[0.5 0; 0 0.5];                   %observer gain

% Receding horizon implementation
x=zeros(dim.nx,T+1);
u_rec=zeros(dim.nu,T);
dhat=zeros(dim.nd,T+1);

%% Simulation 1 MPC feedback control without disturbance
x(:,1)=LTI.x0;
u = zeros(dim.nu, simulation_steps);

for k=1:simulation_steps
    
    %Compute estimated optimal ss
    eqconstraints=eqconstraintsgen(LTI,dim);
    [xr,ur]=optimalss(LTI,dim,weight,[],eqconstraints); 
    
    x_0=x(:,k);
     
    %Solve optimization problem    
    uostar = sdpvar(dim.nu*dim.N,1);                               %define optimization variable
    Constraint=[];                                                 %define constraints
    Objective = 0.5*uostar'*H*uostar+(h*[x_0; xr; ur])'*uostar;    %define cost function
    optimize(Constraint,Objective);                                %solve the problem
    uostar=value(uostar);      

    % Select the first input only
    u(:,k)=uostar(1:dim.nu);
    x(:,k+1) = LTI.A*x(:,k) + LTI.B*Ts*u(:,k);



    % Compute the state/output evolution
    % x(:,k+1)=LTI.A*x_0 + LTI.B*u_rec(:,k);
    % y=LTI.C*x(:,k+1)+LTI.D*u_rec; %LTI.Cd*LTI.d+0.01*randn(dim.ny,1);
    % clear u_uncon
    
    % Update disturbance estimation
    % dhat(:,k+1)=dhat(:,k)+L*(y-LTI.C*x(:,k+1)-dhat(:,k));
 
end

%% Plot simulation 1 results

%PLOTTER(x,t_span);
figure(1)
hold on
plot(t_span(1:end-1), u(1,:));
legend('$\delta$');
title('Control Inputs');

figure(2)
hold on 
plot(t_span(1:end-1), x(1,1:end-1));
%plot(t_span(1:end-1), x(2,1:end-1));
%plot(t_span(1:end-1), x(3,1:end-1));
%plot(t_span(1:end-1), x(4,1:end-1));
%legend('y','ay','r','vy')


