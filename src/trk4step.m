%% Problem 1 - Test Script 

close all; clear; clc;

format long;

% Function that computes right hand sides of ODEs for simple harmonic
% oscillator with unit angular frequency. For use in rk4step, rk4, and rk4ad. 
%
% Governing DE:                                 d2y/dt2 = -y
% Canonical first order dependent variables:    y1 = y, y2 = dydt 
% System of Equations:                          dy1dt = y2, dy2dt = -y1
%
% Inputs
%       t:      Independent variable at current time-step
%       y:      Dependent variables at current time-step (length-n column 
%               vector).
%
% Outputs
%       dydt:  Computes the derivatives of y1 and y2 at the current 
%              time-step (length-n column vector).
function dydt = fcn_sho(t, y)
    dydt = zeros(2,1);
    dydt(1) = y(2);
    dydt(2) = -y(1); 
end

% Function parameters for exact solution of y(t) = sin(t)
x0 = 0; v0 = 1; y0 = [x0; v0];      % Initial conditions 
t0 = 0;                             % Initial time
% Vector of linearly increasing time-step lengths
dt = linspace(0.01, 0.3, 1000);

% Run Runge-Kutta step at various time steps 
yout = zeros(2, length(dt));
for i = 1:length(dt)
    yout(:,i) = rk4step(@fcn_sho, t0, dt(i), y0);
end 

% Calculate the error at each time step length using the known exact solution
errors = abs(yout(1,:) - sin(dt));

% Plot error as a function of dt and compare to C*t^5
hold on;
plot(dt, errors, "Color", 'r');
C = 8.3e-3;
plot(dt, C*dt.^5, "Color", 'b');
title("Magnitude of error vs. time step length dt shown to scale as dt^5");
xlabel("Time step length dt");
ylabel("Magnitude of error");
legend(["Error", "C * t^5"]);
ax = gca; 
ax.FontSize = 12;