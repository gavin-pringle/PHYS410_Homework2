%% Problem 2 - Test Script - Simple Harmonic Oscillator

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
x0 = [0; 1];            % Initial conditions 
t0 = 0; tf = 3*pi;      % Start and end times 

% Vector of output times for each discretization level
tspan6 = linspace(t0, tf, 2^6 + 1);
tspan7 = linspace(t0, tf, 2^7 + 1);
tspan8 = linspace(t0, tf, 2^8 + 1);

% Compute ODE numerical solution at each discretization level
[tout6 xout6] = rk4(@fcn_sho, tspan6, x0);
[tout7 xout7] = rk4(@fcn_sho, tspan7, x0);
[tout8 xout8] = rk4(@fcn_sho, tspan8, x0);

% Plot the solutions at each discretization level 
fig1 = figure(1);
hold on
plot(tout6, xout6(:,1), "LineWidth", 2);
plot(tout7, xout7(:,1), "LineWidth", 2);
plot(tout8, xout8(:,1), "LineWidth", 2);
title("Numerical solutions to SHO ODE at various discretization levels");
xlabel("Independent Variable t");
ylabel("Dependent Variable x");
legend(["l = 6", "l = 7", "l = 8"], 'location', 'best');
ax = gca; 
ax.FontSize = 12;

% Compute the errors at each time step for each discretization level
errors6 = xout6(:,1) - sin(tout6).';
errors7 = xout7(:,1) - sin(tout7).';
errors8 = xout8(:,1) - sin(tout8).';

% Plot the scaled errors for each discretization level 
fig2 = figure(2);
hold on
plot(tout6, errors6, "LineWidth", 2);
plot(tout7, 2^4*errors7, "LineWidth", 2);
plot(tout8, 4^4*errors8, "LineWidth", 2);
grid on 
title({"Scaled errors of numerical solutions to SHO ODE at ", ...
       "various discretization levels"});
xlabel("Independent Variable t");
ylabel("Scaled error");
legend(["error @ l=6", "2^4 * error @ l=7", "4^4 * error @ l=8"], ...
        'location', 'best');
ax = gca; 
ax.FontSize = 12;
