%% Problem 2 - Test Script - Simple Harmonic Oscillator

close all; clear; clc;

format long;

% Function that computes right hand sides of ODEs for simple harmonic
% oscillator with unit angular frequency. For use in rk4step, rk4, and rk4ad. 
%
% Governing DE:                                 y" = -y
% Canonical first order dependent variables:    y1 = y, y2 = y' 
% System of Equations:                          y1' = y2, y2' = -y1
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
t0 = 0; tf = 3*pi;                  % Start and end times 

% Vector of output times for each discretization level
tspan6 = linspace(t0, tf, 2^6 + 1);
tspan7 = linspace(t0, tf, 2^7 + 1);
tspan8 = linspace(t0, tf, 2^8 + 1);

% Compute ODE numerical solution at each discretization level
[tout6 yout6] = rk4(@fcn_sho, tspan6, y0);
[tout7 yout7] = rk4(@fcn_sho, tspan7, y0);
[tout8 yout8] = rk4(@fcn_sho, tspan8, y0);

% Plot the solutions at each 
fig1 = figure(1);
hold on
plot(tout6, yout6(:,1));
plot(tout7, yout7(:,1));
plot(tout8, yout8(:,1));

% Compute the errors at each time step for each discretization level
errors6 = yout6(:,1) - sin(tout6).';
errors7 = yout7(:,1) - sin(tout7).';
errors8 = yout8(:,1) - sin(tout8).';

% Plot the scaled errors for each discretization level 
fig2 = figure(2);
hold on
plot(tout6, errors6);
plot(tout7, 2^4*errors7);
plot(tout8, 4^4*errors8);
grid on 

