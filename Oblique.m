% Oblique Shock wave conditions downstream
% Input 
% Output 

function[M2,P2P1,T2T1,P02P01, beta_deg] = Oblique(M,theta,g,betaguess)
           
        beta = b(M,theta,g,betaguess);
        beta_deg = beta*180/pi;

        % Solve for Mach number behind shock
        % Start by using Mn as input to normal shock solver
        Mn = M*sin(beta);
    
        [M2n, P2P1, T2T1, P02P01] = Normal(Mn,g);
    
    
        M2 = sqrt(M2n^2 + (M*cos(beta))^2/T2T1); 
end


