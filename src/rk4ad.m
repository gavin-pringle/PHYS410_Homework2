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
            % Iteratively compute yf at repeatedly halved dt sizes
            % until reltol is met or floor is reached
            j = 2;
            while dt/(2^j) > floor % Decrease step size by half each iteration
                yc = yf;
                yf = y0;
                for k = 0:2^j - 1 % Number of steps to get to t0 + dt
                    yf = rk4step(fcn, t0 + k*dt/(2^j), dt/(2^j), yf);
                end 
                if abs((yc - yf)/yf) < reltol
                    yout(i,:) = yf.';
                    break;
                end
                j = j + 1;
            end 
            yout(i,:) = yf.';
        end 
    end 

    % Generate array of output values 
    tout = tspan;
end

% NOTES:
% Consider that tspan points might not be equidistant
% RELATIVE tolerance - divide by y
% 
% Can probably add everything into the while / for loop later
%
% y1of4 = rk4step(fcn, t0 + 0*dt/4, dt/4, y0);
% y2of4 = rk4step(fcn, t0 + 1*dt/4, dt/4, y1of4);
% y3of4 = rk4step(fcn, t0 + 2*dt/4, dt/4, y2of4);
% y4of4 = rk4step(fcn, t0 + 3*dt/4, dt/4, y3of4);