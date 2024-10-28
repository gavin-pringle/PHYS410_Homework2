%% Problem 3 - Test Script - Van der Pol oscillator

close all; clear; clc;

format long;

% Function that computes right hand sides of ODEs for Van der Pol
% Oscillator. Following Tsatsos: https://arxiv.org/pdf/0803.1658
%
% Governing DE: x" = -x - a(x^2 - 1)x'
% Canonical first order dependent variables: x1 = x, x2 = x'
% System of Equations: 
%       x1' = x2
%       x2' = -x1 - a(x1^2 - 1)*x2
% 
% Inputs
%       t:      Independent variable at current time-step
%       x:      Dependent variables at current time-step (length-n column 
%               vector).
%
% Outputs
%       dxdt:  Computes the derivatives of x1 and x2 at the current 
%              time-step (length-n column vector).
function dxdt = fcn_vdp(t, x)
    global a;
    dxdt = ones(2,1);
    dxdt(1) = x(2);
    dxdt(2) = -x(1) - a*(x(1)^2 - 1)*x(2);
end

% Function parameters 
x0 = [1; -6];                       % Initial conditions 
tspan = linspace(0.0, 100, 4097);   % Vector of output times
global a; a = 5;                    % Adjustable parameter
reltol = 1.0e-10;                   % Relative tolerance 

% Compute ODE numerical solution 
[tout xout] = rk4ad(@fcn_vdp, tspan, reltol, x0);

% Plot position vs time
fig1 = figure(1);
plot(tout, xout(:,1), "LineWidth", 2, "Color", "#D95319")
title({"Numerical solution of Van der Pol oscillator ODE using rk4ad", ...
       "Position x vs. Time t, Relative tolerance = 1.0e-10"});
xlabel("Independent Variable - Time t");
ylabel("Dependent Variable - Position x");
ax = gca; 
ax.FontSize = 12;

% Plot phase space evolution 
fig2 = figure(2);
plot(xout(:,1), xout(:,2), "LineWidth", 2, "Color", "#D95319")
title({"Phase space evolution of Van der Pol oscillator ODE using rk4ad", ...
       "Velocity dx/dt vs. Position x, Relative tolerance = 1.0e-10"});
xlabel("Position x");
ylabel("Velocity dx/dt");
ax = gca; 
ax.FontSize = 12;
