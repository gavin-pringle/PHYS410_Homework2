%% Problem 3 - Test Script - Simple Harmonic Oscillator

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
x0 = [0; 1];                            % Initial conditions
tspan = linspace(0.0, 3.0 * pi, 65);    % Vector of output times

% Compute ODE numerical solution at each relative tolerance
[tout5 xout5] = rk4ad(@fcn_sho, tspan, 1.0e-5, x0);
[tout7 xout7] = rk4ad(@fcn_sho, tspan, 1.0e-7, x0);
[tout9 xout9] = rk4ad(@fcn_sho, tspan, 1.0e-9, x0);
[tout11 xout11] = rk4ad(@fcn_sho, tspan, 1.0e-11, x0);

% Plot the solutions at each relative tolerance
fig1 = figure(1);
hold on
plot(tout5, xout5(:,1), "LineWidth", 2);
plot(tout7, xout7(:,1), "LineWidth", 2);
plot(tout9, xout9(:,1), "LineWidth", 2);
plot(tout11, xout11(:,1), "LineWidth", 2);
title("Numerical solutions to SHO ODE from rk4ad at various relative tolerances");
xlabel("Independent Variable t");
ylabel("Dependent Variable x");
legend(["reltol = 1.0e-5", "reltol = 1.0e-7", "reltol = 1.0e-9", ...
        "reltol = 1.0e-11"], 'location', 'best');
ax = gca; 
ax.FontSize = 12;

% Compute the errors at each time step for each discretization level
errors5 = xout5(:,1) - sin(tout5).';
errors7 = xout7(:,1) - sin(tout7).';
errors9 = xout9(:,1) - sin(tout9).';
errors11 = xout11(:,1) - sin(tout11).';

% Plot the errors for each relative tolerance
fig2 = figure(2);
hold on
plot(tout5, errors5, "LineWidth", 2);
plot(tout7, errors7, "LineWidth", 2);
plot(tout9, errors9, "LineWidth", 2);
plot(tout11, errors11, "LineWidth", 2);
grid on 
title({"Errors of numerical solutions to SHO ODE at ", ...
       "various rk4ad relative tolerances"});
xlabel("Independent Variable t");
ylabel("Error");
legend(["reltol = 1.0e-5", "reltol = 1.0e-7", "reltol = 1.0e-9", ...
        "reltol = 1.0e-11"], 'location', 'best');
ax = gca; 
ax.FontSize = 12;
