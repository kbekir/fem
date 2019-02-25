function [eps_con,stress_con] = conman(f_c,latrebdet,longrebdet,secdet,s)
% conman is a function for Mander's confined concrete model   
% INPUTS  : f_c(compressive strength of the concrete)
%         : latrebdet(details of the lateral reinforcement)[# of bars,diameter,fyw]
%         : longrebdet(details of the longitudinal reinforcement)
%         : [# of bars,diameter,fyk,eps_su]
%         : secdet(details of the section)[width,length]
%         : s(spacing between the lateral reinforcement bars)
% OUTPUTS : eps_con(strain)
%         : stress_con(stress)


    eps_con(1) = 0;
    stress_con(1) =0; 
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
    counter = 1;

    for j=0.0001:0.00001:eps_cu
        eps_con(counter+1,1) = j;
        r = (E/(E-(f_cc/eps_cc)));
        x = eps_con(counter+1,1)/eps_cc;
        stress_con(counter+1,1) = (f_cc*x*r)/(r-1+x^r);
        counter = counter+1;
    end
end

