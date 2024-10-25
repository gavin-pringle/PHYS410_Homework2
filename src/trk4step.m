%% Problem 1 - Test Script 

close all; clear; clc;

format long;

% Function that computes right hand sides of ODEs for simple harmonic
% oscillator with unit angular frequency. For use in rk4step, rk4, and rk4ad. 
%
% Governing DE:                                 x" = -x
% Canonical first order dependent variables:    x1 = x, x2 = x' 
% System of Equations:                          x1' = x2, x2' = -x1
%
% Inputs
%       t:      Independent variable at current time-step
%       x:      Dependent variables at current time-step (length-n column 
%               vector).
%
% Outputs
%       dxdt:  Computes the derivatives of x1 and x2 at the current 
%              time-step (length-n column vector).
function dxdt = fcn_sho(t, x)
    dxdt = zeros(2,1);
    dxdt(1) = x(2);
    dxdt(2) = -x(1); 
end

% Function parameters for exact solution of x(t) = sin(t)
x0 = [0; 1];    % Initial conditions 
t0 = 0;         % Initial time
% Vector of linearly increasing time-step lengths
dt = linspace(0.01, 0.3, 1000);

% Run Runge-Kutta step at various time steps 
xout = zeros(2, length(dt));
for i = 1:length(dt)
    xout(:,i) = rk4step(@fcn_sho, t0, dt(i), x0);
end 

% Calculate the error at each time step length using the known exact solution
errors = abs(xout(1,:) - sin(dt));

% Plot error as a function of dt and compare to C*t^5
hold on;
plot(dt, errors, "Color", 'r', "LineWidth", 2);
C = 8.3e-3;
plot(dt, C*dt.^5, "--", "Color", 'b', "LineWidth", 2);
title("Magnitude of error vs. time step length dt shown to scale as dt^5");
xlabel("Time step length dt");
ylabel("Magnitude of error");
legend(["Error", "C * t^5"], 'location', 'best');
ax = gca; 
ax.FontSize = 12;