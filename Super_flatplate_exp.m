% Supersonic flat plate pressures
clear; clc; close all;

M = 1.5;
P1 = 81.4;      %[kPa]
T1 = 255.6;     %[K]

theta = 20;     %[°]

g = 1.4;        % ratio of specific heats
M2guess = 2;      % initial guess

[M2,P2P1, T2T1] = PM(theta,g,M,M2guess)


