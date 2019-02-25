function [stress_con] = hysconmander(f_c,latrebdet,longrebdet,secdet,s,eps_hist)
% hysconmander is a function for Mander's confined concrete model for hysteric strains   
% INPUTS  : f_c(compressive strength of the concrete)
%         : latrebdet(details of the lateral reinforcement)[# of bars,diameter,fyw]
%         : longrebdet(details of the longitudinal reinforcement)
%         : [# of bars,diameter,fyk,eps_su]
%         : secdet(details of the section)[width,length]
%         : s(spacing between the lateral reinforcement bars)
%         : eps_hist(given strain history)
% OUTPUTS : stress_con(stress)

    E = 5000*sqrt(f_c);
    eps_c = 0.002;
    fyw = latrebdet(1,3);
    eps_su = longrebdet(1,4);
    area_s = longrebdet(1,1) * (pi()*longrebdet(1,2)^2)/4;
    c = 30; %Concrete cover for all section(mm)
    b_o = secdet(1,1) - (2*c);
    h_o = secdet(1,2) - (2*c);
    area_sx = 2 * b_o * latrebdet(1,2) * latrebdet(1,1);
    area_sy = 2 * h_o * latrebdet(1,2) * latrebdet(1,1);
    ro_x = area_sx/(s*b_o*h_o);
    ro_y = area_sy/(s*h_o*b_o);
    ro_s = ro_x+ro_y;

    k_e = (1-((2*h_o+2*b_o)/(6*b_o*h_o)))*(1-(s/(2*b_o)))*(1-(s/(2*h_o)))*(1-((area_s)/(b_o*h_o)))^-1;
    fe_x = k_e*ro_x*fyw;
    fe_y = k_e*ro_y*fyw;
    fe = (fe_x+fe_y)/2;
    lambda_c = 2.254*sqrt(1+(7.94*(fe/f_c))) - 2*(fe/f_c) -1.254;

    f_cc = lambda_c*f_c;
    eps_cc = eps_c*(1+5*(lambda_c-1));
    eps_cu = 0.004 + ((1.4*ro_s*fyw*eps_su)/f_cc);

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
            stress_cur = (f_cc*x*r)/(r-1+x^r);
        end
        
        if eps_r > eps_cur
            if eps_cur-eps_r == 0
                stress_con(i,1) = 0;
            elseif eps_cur-eps_r ~= 0
                stress_con(i,1) = stress_cur * (1+((eps_cur-eps_r)/(eps_r-eps_p)));
            end
        elseif eps_r <= eps_cur
            stress_con(i,1) = stress_cur;
        end
        
        if eps_cur < 0
            stress_con(i,1) = 0;
        end
        
        if stress_con(i,1) < 0
            stress_con(i,1) = 0;
        end
    end
end
