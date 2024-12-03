% Beta angle for Oblique solver
% Input incoming Mach, theta [deg], ratio of specific heats, betaguess
% [deg]
% Output Beta [rad]

function beta = b(M,theta,g,betaguess)
    beta = pi/180 * betaguess;
    theta = pi/180 * theta;
    diff = 100;
    
    while abs(diff) > 1E-7
        diff = 2*(M^2 * sin(beta)^2 -1)/(tan(beta)*(2+M^2*(g + cos(2*beta)))) - tan(theta);
        beta = beta - .1*diff;
    end
    
end
       

    