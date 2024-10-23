%% Problem 1 - Single Fourth Order Runge-Kutta Step

% Inputs
%       fcn: Function handle for right hand sides of ODEs (returns
%       length-n column vector).
%       t0: Initial value of independent variable.
%       dt: Time step.
%       y0: Initial values (length-n column vector).
%
% Output
%       yout: Final values (length-n column vector)
function yout = rk4step(fcn, t0, dt, y0)
    yout = 0;