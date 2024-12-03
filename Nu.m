%Prandtl angle given Mach and gamma. For expansion shock waves.
% Inputs : Mach number (M) , ratio of specific heats (g)
% Output: angle in radians

function[r] = Nu(M,g)
    r = sqrt((g+1)/(g-1)) * atan(sqrt(((g-1)/(g+1))*(M^2-1))) - atan(sqrt((M^2 - 1)));
end
