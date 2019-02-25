clear all
clc

% Test code for confined and unconfined mander concrete model

% Unconfined mander concrete model needs only the compressive strength of
% concrete(f_c) and the strain history as input

loading_type = 1; 						%Give the loading type details (1 for backbone, 2 for hysteric)

f_c = 30;

if loading_type == 1
    [eps_uncon,stress_uncon] = unconman(f_c);
    plot(eps_uncon,stress_uncon)
    grid on;
    hold on;
    
    % Confined mander concrete model needs the compressive strength of
    % concrete(f_c), the strain history(eps_hist), details for lateral reinforcement,
    % details for longitudinal reinforcement, spacing between the lateral
    % reinforcement bars and the details of the section as input

    latrebdet = [1,12,420]; 			%Details for lateral reinforcement [# of bar,diameter(mm),fyw(mpa)]
    longrebdet = [4,18,420,0.12];   	%Details for longitudinal reinforcement [# of bar,diameter(mm),fyk(mpa),eps_su]
    secdet = [300,500];  				%Details for section [width(mm),length(mm)]
    s = 50; 							%Spacing for lateral reinforcement(mm)
    [eps_con,stress_con] = conman(f_c,latrebdet,longrebdet,secdet,s);
    
    plot(eps_con,stress_con)
    grid on;

elseif loading_type == 2

	%In this section of the code, the program reads the given strain history
	%from the given file. Then, locate the strain data to the eps_hist matrix

    fileID = fopen('epshist.txt','r'); 	%The name of the file, which contains the strain data, can be written 
    formatSpec = '%f'; 					%Format of the strain data
    eps_hist = fscanf(fileID,formatSpec);

    [stress_uncon] = hysunconmander(f_c,eps_hist);

    plot(eps_hist,stress_uncon)
    grid on;
    hold on;

    % Confined mander concrete model needs the compressive strength of
    % concrete(f_c), the strain history(eps_hist), details for lateral reinforcement,
    % details for longitudinal reinforcement, spacing between the lateral
    % reinforcement bars and the details of the section as input

    latrebdet = [1,12,420]; 			%Details for lateral reinforcement [# of bar,diameter(mm),fyw(mpa)]
    longrebdet = [4,18,420,0.12]; 		%Details for longitudinal reinforcement [# of bar,diameter(mm),fyk(mpa),eps_su]
    secdet = [300,500]; 				%Details for section [width(mm),length(mm)]
    s = 50; 							%Spacing for lateral reinforcement(mm)

    [stress_con] = hysconmander(f_c,latrebdet,longrebdet,secdet,s,eps_hist);
    plot(eps_hist,stress_con)
    grid on;
    xlim ([-0.5*10^-3 6*10^-3]);
    ylim ([-5 35]);

end
