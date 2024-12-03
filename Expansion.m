%Mach number downstream of Prandtl-Meyer Expansion Fan
% Input theta in degrees
% Output M2, ratio of P2/P1, T2/T1

function[M2, Pratio, Tratio] = Expansion(theta,g,M1,M2guess)
    theta = theta * pi/180;
    nu1 = Nu(M1,g);
    nu2 = theta + nu1;
    for j = 1:100
        M2 = (theta +  nu1- Nu(M2guess,g))/dNu(M2guess,g) + M2guess;
        M2guess = M2;
    end

    Tratio = (1 + 0.5*(g-1)*M1^2)/(1 + 0.5*(g-1)*M2^2);
    Pratio = Tratio^(g/(g-1));
end
