%% Problem 3 - Adaptive Fourth Order Runge-Kutta System of ODEs Integrator

%
%
% Inputs
%       fcn:    Function handle for right hand sides of ODEs (returns
%               length-n column vector)
%       tspan:  Vector of output times (length nout vector).
%       reltol: Relative tolerance parameter.
%       y0:     Initial values (length-n column vector).
%
% Outputs
%       tout:   Output times (length-nout column vector, elements
%               identical to tspan).
%       yout:   Output values (nout x n array. The ith column of yout
%               contains the nout values of the ith dependent variable).
function [tout yout] = rk4ad(fcn, tspan, reltol, y0)
    % Number of equations in ODE system
    n = max(size(y0));
    % Number of time-steps
    nout = max(size(tspan));

    % Lower bound on step size
    floor = 1.0e-4;

    % Initialize array for output values 
    yout = zeros(nout, n);
    yout(1,:) = y0.';

    % Integrate ODE 
    for i = 2:nout
        yout(i,:) = rk4step(fcn, tspan(i-1), dt, yout(i-1,:).').';
    end 
end

% NOTES:
% Consider that tspan points might not be equidistant
% tout may have more entries than tspan 
% 