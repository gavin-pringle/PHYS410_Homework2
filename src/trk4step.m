%% Problem 1 - Test Script 

close all; clear; clc;

format long

% Function handle for right hand sides of ODEs (returns length-n column vector).
function dydt = fcn_sho(t, y)
    % fcn_sho: Derivatives function for simple harmonic motion with unit
    % angular frequency. Note that the function must return a COLUMN
    % vector.

    % d2y/dt2 = -y
    % y1 = y, y2 = dy1dt, dy2dt = -y1
    
    % The exact solution to this ODE is y(t) = sin(t)

    dydt = zeros(2,1);
    dydt(1) = y(2);
    dydt(2) = -y(1); 
end

% Function parameters 
x0 = 0; v0 = 1; y0 = [x0; v0];      % Initial conditions 
t0 = 0;                             % Initial time
dt = linspace(0.01, 0.3, 1000);   % Time-steps

% Run Runge-Kutta step at various time steps 
yout = zeros(2, length(dt));
for i = 1:length(dt)
    yout(:,i) = rk4step(@fcn_sho, t0, dt(i), y0);
end 

% Calculate the error at each time step length using the known exact solution
errors = abs(yout(1,:) - sin(dt));

plot(dt, errors, "LineWidth", 1, "Color", 'r');
hold on 
a = 0.83e-2;
plot(dt,a*dt.^5, "LineWidth", 1, "Color", 'b')
