%% Problem 2 - Runge-Kutta 

% Inputs
%       fcn: Function handle for right hand sides of ODEs (returns
%       length-n column vector)
%       tspan: Vector of output times (length nout).
%       y0: Initial values (length-n column vector).
%
% Outputs
%       tout: Vector of output times.
%       yout: Output values (nout x n array. The ith column of yout
%       contains the nout values of the ith dependent variable).
function [tout yout] = rk4(fcn, tspan, y0)
    tout = 0;
    yout = 0;