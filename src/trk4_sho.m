%% Problem 2 - Test Script - Simple Harmonic Oscillator

close all; clear; clc;

format long;

% Function that computes right hand sides of ODEs for simple harmonic
% oscillator with unit angular frequency. For use in rk4step, rk4, and rk4ad. 
%
% Governing DE:                                 d2y/dt2 = -y
% Canonical first order dependent variables:    y1 = y, y2 = dydt 
% System of Equations:                          dy1dt = y2, dy2dt = -y1
%
% Inputs
%       t:      Independent variable at current time-step
%       y:      Dependent variables at current time-step (length-n column 
%               vector).
%
% Outputs
%       dydt:  Computes the derivatives of y1 and y2 at the current 
%              time-step (length-n column vector).
function dydt = fcn_sho(t, y)
    dydt = zeros(2,1);
    dydt(1) = y(2);
    dydt(2) = -y(1); 
end

% Function parameters for exact solution of y(t) = sin(t)
x0 = 0; v0 = 1; y0 = [x0; v0];      % Initial conditions 
t0 = 0;                             % Initial time