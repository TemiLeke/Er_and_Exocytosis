function [Jpmca, Jserca, Jvgcc, Jdiff_vgcc, Jipr, Jdiff,...
           Jleak, Jin, JCa_vgcc_Ca_ipr] = Fluxes(Cac, Ca_vgcc, Cer, Ca_ipr, IP3, ...
                                           IcaPQ, Po_ip3r, coupling_condition, cell_condition)

global N_vgcc

%% PMCA FLUX
    
    Vpmca = 3.195;   % uM/ms   Maximum capacity of plasma pump
    Kpmca = 0.5;     % uM      Half-maximal activating Cac of plasma pump     
    nP = 2.0;        %         Hill coefficient of plasma pump

    Jpmca = Vpmca*(Cac^nP)/(Kpmca^nP + Cac^nP);           % uM/ms  Flux through plasma pump
    
%% SERCA Flux 

    Vs = 10.0;       % uM/ms   Maximum capacity of SERCA
    Ks = 0.26;       % uM      SERCA half-maximal activating Cac
    ns = 1.75;       %         Hill coefficient of SERCA
    
    Jserca = Vs * (Cac^ns) / (Ks^ns + Cac^ns);            % uM/ms  Uptake of Ca into the ER by SERCA pumps

%% VGCC Flux
    
%   Parameters (some modified) from Schikorski. et.
%   al. 1997, Holderith et. al. 2012, and Ermolyuk. et. al 2013, 

    z = 2;                                      %         valence of Ca ion
    F = 96485.33;                               % C/mole  Faraday's constant, unit: coul/mole
    Area = 3.8489e-09;                          % cm^2    Bouton Membrane area in cm^2
    Vol_tmnal = 1.22e-16;                       % Litre   Presynaptic terminal Volume (Assuming a spherical Bouton shape)
    Kdiff_vgcc = 0.071;                         % /ms     Rate of diffusion from VGCC nanodomain
    cluster_radius = 25e-03;                    % um
    cluster_area = (pi)*(cluster_radius)^2;     % um^2
    active_zone_area = 0.04;                    % um^2
    active_zone_number = 1.3;                   %         Number of active zones
    
    channel_density = N_vgcc / (active_zone_number * active_zone_area);
    
    Ica = channel_density * cluster_area * IcaPQ;  
    
    % Macroscopic current density through an ensemble of open channels where channel_density is crucial:
    Jvgcc = (-Ica/(z*F*Vol_tmnal))*1e-03;         % uM/ms   Calcium influx from the extracellular space to cytosol(with Mitoc 2.03)
    Jdiff_vgcc = Kdiff_vgcc * (Ca_vgcc - Cac);    % uM/ms   Diffusion from VGCC cluster nanodomain to the cytoplasm
    
%% Receptors -- IP3 and RYR

    KIPR = 5;                                      % /ms  IP3R flux coefficient 
    kdiff = 10;                                    % /ms  Ca diffusional flux coefficient
    
    Jipr  = KIPR * Po_ip3r * (Cer - Ca_ipr);       % uM/ms  Flux through the IP3R
    Jdiff = kdiff * (Ca_ipr - Cac);                % uM/ms  Diffusion from ER cluster nanodomain to the cytoplasm
    
 %% Leak Fluxes -- ER and Plasma Membrane

    kleak = 0.0022;                   % |   /ms     |	ER leak flux coefficient
    Jleakin = 0.03115;                % |   uM/ms   |	Plasma membrane leak influx
    Vleakin = 0.2;                    % |   /ms     |	IP3 In Leak flux coefficient

    Jleak = kleak * (Cer - Cac);      % |   uM/ms   |	Baground leak from the ER into the cytoplasm    
    Jin = Jleakin + Vleakin*IP3; 

 %% Flux from VGCC Nano Domain to IP3R Nanodomain
    if coupling_condition == "Higher_Coupling_AD"
        if cell_condition == "WT"
            Kc = 20; 
            Vc = 118;
            K_hat = 5;
        elseif cell_condition == "AD"
            Kc = 10; 
            Vc = 118;
            K_hat = 15;
        end
    elseif coupling_condition == "Higher_Coupling_WT"
        if cell_condition == "WT"
            Kc = 10; 
            Vc = 118;
            K_hat = 15;
        elseif cell_condition == "AD"
            Kc = 20; 
            Vc = 118;
            K_hat = 5;
        end
    elseif coupling_condition == "Same_Coupling"
            Kc = 10; 
            Vc = 118;
            K_hat = 15;
    end
    
    JCa_vgcc_Ca_ipr = Vc * (Ca_vgcc^2 - K_hat*Ca_ipr^2)/(Kc^2 + Ca_vgcc^2);
    