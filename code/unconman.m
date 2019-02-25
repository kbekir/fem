function [eps_uncon,stress_uncon] = unconman(f_c)
% unconman is a function for Mander's unconfined concrete model   
% INPUTS  : f_c(compressive strength of the concrete)
% OUTPUTS : eps_con(strain)
%         : stress_con(stress)


    eps_uncon(1) = 0;
    stress_uncon(1) =0; 
    E = 5000*sqrt(f_c);
    eps_c = 0.002;
    eps_u = 0.005;
    counter = 1;

    for i=0.0001:0.00001:eps_u
        eps_uncon(counter+1,1) = i;
        r = (E/(E-(f_c/eps_c)));
        x = eps_uncon(counter+1,1)/eps_c;

        if eps_uncon(counter+1,1)<=2*eps_c
           stress_uncon(counter+1,1) = ((f_c)*x*r)/(r-1+x.^(r));
        end

        if 2*eps_c<eps_uncon(counter+1,1)
           stress_uncon(counter+1,1) = (((2*f_c*r)/(r-1+2^(r)))*(eps_u-eps_uncon(counter+1,1))/(eps_u-2*eps_c));
        end

        counter = counter+1;
    end
    
end

