%% Problem 2 - Runge-Kutta System of ODEs Integrator

% Function that numerically computes the solution to a system of ODEs
% over a given period of time using a fourth-order Runge-Kutta method.
%
% Inputs
%       fcn:    Function handle for right hand sides of ODEs (returns
%               length-n column vector)
%       tspan:  Vector of output times (length nout).
%       y0:     Initial values (length-n column vector).
%
% Outputs
%       tout:   Vector of output times.
%       yout:   Output values (nout x n array. The ith column of yout
%               contains the nout values of the ith dependent variable).
function [tout yout] = rk4(fcn, tspan, y0)
    % Number of equations in ODE system
    n = max(size(y0));
    % Number of time-steps
    nout = max(size(tspan));
    % Step size
    dt = tspan(2) - tspan(1);

    % Initialize array for output values 
    yout = zeros(nout, n);
    yout(1,:) = y0.';

    % Integrate ODE 
    for i = 2:nout
        yout(i,:) = rk4step(fcn, tspan(i-1), dt, yout(i-1,:).').';
    end 

    % Generate array of output values 
    tout = tspan;
end