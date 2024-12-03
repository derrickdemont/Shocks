% Oblique Shock Conditions
clear; clc; close all;

M1 = 2;
g = 1.4;        % ratio of specific heats
betaguess = 50;      % initial guess
theta = 15;

[M2,P2P1,T2T1,P02P01, beta_deg] = Oblique(M1,theta,g,betaguess);

%% Create figure
x = [0 1 0 cosd(theta) 0 cosd(beta_deg)];
y = [0 0 0 sind(theta) 0 sind(beta_deg)];

Linewidth = 1.5;
f=figure;
plot(x(1:2),y(1:2), x(3:4), y(3:4), 'Color', [0 0.4470 0.7410], 'LineWidth', Linewidth)
hold on
plot(x(5:6), y(5:6), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', Linewidth)
xlim([-0.35 1.1])
ylim([-0.5 1])

% Labels
r = 0.3;
bangs = 0:.01:beta_deg;
xb = r*cosd(bangs);
yb = r*sind(bangs);
plot(xb,yb, 'Color', 'k', 'LineWidth', Linewidth)

r = 0.5;
tangs = 0:.01:theta;
xt = r*cosd(tangs);
yt = r*sind(tangs);
plot(xt,yt, 'Color', 'k', 'LineWidth', Linewidth)


r = 0.7;
text(r*cosd(theta/2),r*sind(theta/2),strcat('θ = ', num2str(theta), '°'))
r = 0.5;
text(r*cosd((beta_deg-theta)*0.6+theta),r*sind((beta_deg-theta)*0.6+theta),strcat('β = ', num2str(beta_deg), '°'))

% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);

dim = [.2 .5 .3 .3];
% str = strcat('M1 = ', num2str(M1), '   ', '  M2 = ', num2str(M2));
% str = str + newline + strcat('P2/P1 = ', num2str(P2P1));
% str = str + newline + strcat('T2/T1 = ', num2str(T2T1));
% str = str + newline + strcat('P02/P01 = ', num2str(P02P01));

str = {strcat('M1 = ', num2str(M1), '   ', '  M2 = ', num2str(M2)),strcat('P2/P1 = ', num2str(P2P1)),strcat('T2/T1 = ', num2str(T2T1)),strcat('P02/P01 = ', num2str(P02P01))};

annotation('textbox',dim,'String',str,'FitBoxToText','on');





%% Functions

function[M2,P2P1,T2T1,P02P01, beta_deg] = Oblique(M,theta,g,betaguess)
           
        beta = b(M,theta,g,betaguess);
        beta_deg = beta*180/pi;

        % Solve for Mach number behind shock
        % Start by using Mn as input to normal shock solver
        Mn = M*sin(beta);
    
        [M2n, P2P1, T2T1, P02P01] = Normal(Mn,g);
    
    
        M2 = sqrt(M2n^2 + (M*cos(beta))^2/T2T1); 
end

function beta = b(M,theta,g,betaguess)
    beta = pi/180 * betaguess;
    theta = pi/180 * theta;
    diff = 100;
    
    while abs(diff) > 1E-7
        diff = 2*(M^2 * sin(beta)^2 -1)/(tan(beta)*(2+M^2*(g + cos(2*beta)))) - tan(theta);
        beta = beta - .1*diff;
    end
    
end
     

function [M2, P2P1, T2T1, P02P01] = Normal(M1,g)
    M2 = sqrt((1+ 0.5*(g-1)*M1^2)/(g*M1^2 - 0.5*(g-1)));
    P2P1 = 1 + 2*g*(M1^2 - 1)/(g+1);
    T2T1 = (1 + 2*g*(M1^2 - 1)/(g+1))*((2 + (g-1)*M1^2)/((g+1)*M1^2));

    A = 2/((g+1)*(g*M1^2 - 0.5*(g-1))^(1/(g-1)));
    B = (0.5*(g+1)*M1)^2;
    C = 1 + 0.5*(g-1)*M1^2;
    P02P01 = A*(B/C)^(g/(g-1));
    
end
       

    


