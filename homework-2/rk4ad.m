%% Problem 3 - Adaptive Fourth Order Runge-Kutta System of ODEs Integrator

% Function that numerically computes the solution to a system of ODEs
% over a given period of time using a fourth-order Runge-Kutta method
% with adaptive steps sizes to ensure a relative tolerance is reached.
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
        tprev = tspan(i-1);
        yprev = yout(i-1,:).';
        dt = tspan(i) - tspan(i-1);

        % Compute fine and coarse approximations for y(tprev + dt)
        yc = rk4step(fcn, tprev, dt, yprev);
        if dt/2 < floor
            % If fine step is lower than floor, cannot narrow down any further
            yout(i,:) = yc.';
            continue;
        end 
        yhalf = rk4step(fcn, tprev, dt/2, yprev);
        yf = rk4step(fcn, tprev + dt/2, dt/2, yhalf);

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
                yf = yprev;
                for k = 0:2^j - 1 % Number of steps to get to tprev + dt
                    yf = rk4step(fcn, tprev + k*dt/(2^j), dt/(2^j), yf);
                end 

                if abs((yc - yf)/yf) < reltol
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
