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
    
end