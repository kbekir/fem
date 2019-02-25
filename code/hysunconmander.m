function [stress_uncon] = hysunconmander(f_c,eps_hist)
% hysconmander is a function for Mander's unconfined concrete model   
% INPUTS  : f_c(compressive strength of the concrete)
%         : eps_hist(given strain history)
% OUTPUTS : stress_con(stress)

    E = 5000*sqrt(f_c);
    eps_c = 0.002;
    eps_u = 0.005;

    for i=1:size(eps_hist,1)
        eps_cur = eps_hist(i,1);
        
        if i==1
            eps_r = eps_cur;
        elseif i>1
            eps_r = max(eps_r,eps_cur);
        end
        
        eps_ratio = eps_r/eps_c;
        
        if eps_ratio < 2
            eps_p = eps_c*(((0.145)*(eps_ratio)^2) + (0.13*eps_ratio));
        elseif eps_ratio >= 2
            eps_p = eps_c*((0.707*(eps_ratio-2) + 0.834));
        end
       
        if eps_cur < 0
            stress_cur = 0;
        elseif eps_cur == 0
            stress_cur = 0;
        elseif eps_cur > 0
            r = (E/(E-(f_c/eps_c)));
            x = eps_cur/eps_c;

            if eps_cur <= 2*eps_c
               stress_cur = ((f_c)*x*r)/(r-1+x.^(r));
            elseif 2*eps_c<eps_cur
               stress_cur = (((2*f_c*r)/(r-1+2^(r)))*(eps_u-eps_cur)/(eps_u-2*eps_c));
            end

        if eps_r > eps_cur
            if eps_cur-eps_r == 0
                stress_uncon(i,1) = 0;
            elseif eps_cur-eps_r ~= 0
                stress_uncon(i,1) = stress_cur * (1+((eps_cur-eps_r)/(eps_r-eps_p)));
            end
        elseif eps_r <= eps_cur
            stress_uncon(i,1) = stress_cur;
        end
        
        if eps_cur < 0
            stress_uncon(i,1) = 0;
        end
        
        if stress_uncon(i,1) < 0
            stress_uncon(i,1) = 0;
        end   
    end    
end
