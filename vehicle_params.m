%% Vehicle model parameters

g = 9.81;           % gravity acceleration m/s^2
l_f = 1.40;         % distance from front axle to CoG, m 
l_r = 1.45;         % distance from rear axle to CoG, m 
L = l_f + l_f;      % wheelbase, m 
m = 1950;           % vehicle mass, kg
Izz = 3500;         % inertia moment of vehicle about vertical axis, kg*m^2
Ca_f = 2*92000;     % front axle cornering stiffness, N/rad 
Ca_r = 2*97000;     % rear axle cornering stiffness, N/rad 
usteer = 17.5;        % steering ratio

