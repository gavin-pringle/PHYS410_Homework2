%% Problem 2 - Test Script - Van der Pol oscillator

close all; clear; clc;

format long;

% Function that computes right hand sides of ODEs for Van der Pol
% Oscillator. Following Tsatsos: https://arxiv.org/pdf/0803.1658
%
% Governing DE: y" = -y - a(y^2 - 1)y' + b*sin(w*t)
% Canonical first order dependent variables: y1 = y, y2 = y'
% System of Equations: 
%       y1' = y2
%       y2' = -y1 - a(y1^2 - 1)*y2 + b*sin(w*t)
% 
% Inputs
%       t:      Independent variable at current time-step
%       y:      Dependent variables at current time-step (length-n column 
%               vector).
%
% Outputs
%       dydt:  Computes the derivatives of y1 and y2 at the current 
%              time-step (length-n column vector).
function dydt = fcn_vdp(t, y)
    global a b omega;
    dydt = ones(2,1);
    dydt(1) = y(2);
    dydt(2) = -y(1) - a*(y(1)^2 - 1)*y(2) + b*sin(omega*t);
end

