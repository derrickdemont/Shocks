%Derivative of Prandtl-Meyer angle Nu(M) with respect to M. For expansion shock waves
% Inputs : Mach number (M) , ratio of specific heats (g)
% Output : deriv of nu wrt M

function[r] = dNu(M,g)
    r = (1/M)*(sqrt(M^2 - 1)) / (1+(M^2)*((g-1)/2));
end
