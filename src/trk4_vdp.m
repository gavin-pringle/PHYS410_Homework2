%% Problem 2 - Test Script - Van der Pol oscillator

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
x0 = [0; -6];       % Initial conditions 
t0 = 0; tf = 100;   % Start and end times 
a = 5;              % Adjustable parameter

% Discretization level 
level = 12; 
tspan = linspace(t0, tf, 2^level + 1);

% Compute ODE numerical solution 
[tout xout] = rk4(@fcn_vdp, tspan, x0);

% Plot position vs time
fig1 = figure(1);
plot(tout, xout, "LineWidth", 2)
title("Numerical solution of Van der Pol oscillator ODE - Position x vs. Time t");
xlabel("Independent Variable - Time t");
ylabel("Dependent Variable - Position x");
ax = gca; 
ax.FontSize = 12;

% Plot phase space evolution 

