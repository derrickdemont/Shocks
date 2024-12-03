% Supersonic flat plate pressures
clear; clc; close all;

M = 3;
Pabs = 650;
AOA = 20;
P0 = 0*AOA;
P = P0;

g = 1.4;        % ratio of specific heats
betaguess = 55;      % initial guess

for i=1:length(AOA)

        theta = AOA(i);     % angle of attack
        
        
        beta = b(M,theta,g,betaguess);
        beta_deg = beta*180/pi;

        % Solve for Mach number behind shock
        % Start by using Mn as input to normal shock solver
        Mn = M*sin(beta);
    
        [M2n, P2P1, T2T1, P02P01] = Normal(Mn,g);
    
    
        M2 = sqrt(M2n^2 + (M*cos(beta))^2/T2T1);
        
        
 
        Pratio = 1 + ((2*g)/(g+1))*((M*sin(beta))^2 - 1);
        P(i) = Pabs *Pratio       % Impinging side pressure
        

        % Dynamic pressure
        q(i) = 0.5*P(i)*g*M2^2
        
        % Calculate stagnation pressure in crack given local pressure, Mach
        
        P0(i) = P(i)*(1 + 0.5*(g-1)*M2^2)^(g/(g-1))     % [psi]

end

P0_inf = Pabs*(1 + 0.5*(g-1)*M^2)^(g/(g-1));
a = 5*pi/180;
M2NASA = sqrt(((g-1)*M^2*sin(beta)^2 + 2)/(sin(beta-a)^2 * (2*g*M^2*sin(beta)^2 - (g-1))))


