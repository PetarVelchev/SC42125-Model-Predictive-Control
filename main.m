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
% state vector x = [v r]
% control vector u = [delta]
% output vector y = [ay r beta v]
% delta - road steering angle, rad
% ay - lateral acceleration, m/s^2
% r - yaw rate, rad/s
% v - lateral velocity, m/s
% beta - sideslip angle, rad

A_ss = -[Ct / (m * Vx), Cs / (m * Vx) + Vx; Cs / (Izz * Vx), Cq2 / (Izz * Vx)];
B_ss = [Ca_f / m; l_f * Ca_f / Izz];
C_ss = [-Ct / (m * Vx), -Cs / (m * Vx); 0, 1; 1/Vx, 0; 1 0];
D_ss = [Ca_f / m; 0; 0; 0];

% Define transfer function
sys_cont = ss(A_ss,B_ss,C_ss,D_ss);
systf_cont = tf(sys_cont);

%% Discrete system dynamics 
T = 10; %Simulation Time
Ts = 0.1; %Sampling Time
t_span = 0:Ts:T; %span vector
x0 = [Vx, 0]';
simulation_steps = T/Ts; 

sys_disc = c2d(sys_cont, Ts);