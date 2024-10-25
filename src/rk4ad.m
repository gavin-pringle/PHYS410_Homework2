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
        % Compute coarse rk4step arguments  
        t0 = tspan(i-1);
        y0 = yout(i-1,:).';
        dt = tspan(i) - tspan(i-1);

        % Compute fine and coarse appromations for y(t0 + dt)
        yc = rk4step(fcn, t0, dt, y0);
        if dt/2 < floor
            % If fine step is lower than floor, cannot narrow down any further
            yout(i,:) = yc.';
            continue;
        end 
        yhalf = rk4step(fcn, t0, dt/2, y0);
        yf = rk4step(fcn, t0 + dt/2, dt/2, yhalf);

        % Check if error meets relative tolerance parameter
        if abs((yc - yf)/yf) < reltol
            yout(i,:) = yf.';
            continue;
        else 
            % Iteratively halve dt until reltol is met or floor is reached
        end 
    end 
end

% NOTES:
% Consider that tspan points might not be equidistant
% RELATIVE tolerance - divide by y
% 