function p = parameters_homogeneous(varargin)


    if nargin==0
        p.R = 8.314; % Gas constant (KJ/Kmol K)
        p.g = 9.81; % Accelaration due to gravity (m/s2)
        p.ST = 22.5e-3; % Surface tesnison (N/m)

        %Temperature
        p.T_ref = 50 +273.15; % reference Temperature (K) 
        p.T = 55 + 273.15; % Reactor Temeprature (K)

        %Pressure
        p.A_MeOH = 5.2041; % Antoine constant A for methanol
        p.B_MeOH = 1581.341; % Antoine constant B for Methanol
        p.C_MeOH = -33.5; % Antoines constant C for methanol 
        p.vap_press_MeOH = 10^(p.A_MeOH -(p.B_MeOH/(p.T + p.C_MeOH))); % vapor Pressure of methanol at 65C (bar)
        p.P_R = 1.5*p.vap_press_MeOH; %Sytem presure 50% higher than vapor pressure of Methanol (N/m2) taken 

        %% TG PHASE/Continuous Phase
        p.rho_TG = 918.8; % Density of TG, (kg/m3)
        p.MW_TG = 794.326; % Molecular weight of TG (g/mol)

        %flow rates of TG
        p.mdot_TG = 1.4/320; % Mass Flow rate of TG into reactor (kg/s)
        p.F_TG = (p.mdot_TG*1000) / p.MW_TG; % Molar Flow rate of TG phase (mol/s);
        p.C0_TG = (p.rho_TG/p.MW_TG)*1000; % Inital concentration of TG  in the reactor (mol/m3_TG)
        p.Fv_TG = p.F_TG/p.C0_TG; % Volumetric flow rate of TG phase(m3_TG/s)

        %viscosity
        p.mu_TG_ref = exp(-13.2722 + (2792.919/p.T_ref)); % Dynamic viscosity of TG phase at reference temperature (Pa.s)
        p.mu_TG = exp(-13.2722 + (2792.919/p.T)); % Dynamic viscosity of TG phase at reactor temperature (Pa.s)
        p.mu_k_TG = p.mu_TG/p.rho_TG; % Kinematic viscosity of TG phase (m2/s)

        %diffusivity
        p.D_TG_ref = [6.1e-11, 5.8e-11, 5.7e-11, 5.6e-11, 3.7e-11, 7.8e-10]; % Reference Diffusivities for TG phase componenets(m2/s)
        p.D_TG = p.D_TG_ref*(p.T/p.T_ref)*(p.mu_TG_ref/p.mu_TG); % Diffusivities for TG phase componenets at specified temperature (m2/s)

        %% METHANOL PHASE/Disperesed Phase
        p.MW_MeOH =32; % Molecular weight of Methanol (g/mol)
        p.MW_Na = 40; % Molecular weight of NaOH (kg/kmol)
        p.rho_MeOH = 792; % Density of dispersed phase (kg/m3)
        p.ratio = 6; % Molar Ratio of Methanol to TG

        %flow rate
        p.F_MeOH = p.ratio*p.F_TG; % Molar flow of methanol phase(mol/s)
        p.mdot_MeOH = (p.MW_MeOH*p.F_MeOH)/1000; % Mass flow rate of Methanol (kg/s)
        p.C0_MeOH = (p.rho_MeOH/p.MW_MeOH)*1000; % Inital concentration of  methanol in Methanol phase (mol/m3_ME)
        p.Fv_MeOH = p.F_MeOH/p.C0_MeOH; % Volumetric flow rate of methanol phase (m3_ME/s)

        %viscosity
        p.mu_MeOH_ref = exp(-11.916 + (1315.706/p.T_ref)); % Dynamic viscosity of Methanol phase at reference temperature (Pa.s)
        p.mu_MeOH = exp(-11.916 + (1315.706/p.T)); % Dynamic viscosity of Methanol phase at reactor temperature (Pa.s)
        p.mu_k_MeOH = p.mu_MeOH/p.rho_MeOH; % Kinematic viscosity of Methanol phase(m2/s)
        
        %Diffusivity
        p.D_ME_ref = [9.2e-11, 8.4e-11, 7.4e-11, 6.2e-11, 4.6e-11, 1.6e-9]; % Reference Diffusivities for TG phase componenets(m2/s)
        p.D_ME = p.D_ME_ref*(p.T/p.T_ref)*(p.mu_MeOH_ref/p.mu_MeOH); % Diffusivities for TG phase componenets at specified temperature (m2/s)


        %%Homogeneous system

        % DISTRIBUTION COEFFICIENT
        p.m_DC = [150, 105, 54, 0.0062, 47, 0.0031]; % m3_MeOH/m3_TG (m_TG,m_DG,m_MG,m_G,m_FAME,m_MeOH)


        % REACTOR and REACTION PARAMETERS
        p.h = 2e-3; % Gap distance in RS-SDR (m)
        p.rd = 65e-3; % spinner dia (m)
        p.V_r = 6.2e-5; % Volume of reactor(where reaction takes place) (m3_R)
        p.Fv_Tot = p.Fv_TG + p.Fv_MeOH; % Total Volumetric Flow rate; (m3_R/s)
        p.eps_TG = p.Fv_TG /p.Fv_Tot; % Hold up of contiouns pahse in reactor (m3_TG/m3_R)
        p.eps_MeOH = p.Fv_MeOH /p.Fv_Tot; % Hold up of Methanol pahse in reactor (m3_MeOH/m3_R)
        p.V_TG = p.eps_TG * p.V_r; % Volume of TG Phase (m3_TG)
        p.V_MeOH = p.eps_MeOH * p.V_r; % Volume of Methanol Phase (m3_ME)

        % REACTION KINETICS FOR HOMOGENEOUS REACTIONS
        p.k_hom_ref = [0.1172/60, 0.1680/60, 0.1372/60, 0.0565/60, 0.0870/60, 0.018/60]; % Rate constants at 50 C for all 6 reactions in Homogeneous case (m3/mol.s)
        p.EA_hom = [85.9e3, 82.8e3, 17.3e3, 90.3e3, 38.1e3, 13.5e3]; % Activation energy (J/mol)
        p.k_hom = p.k_hom_ref .* exp((-p.EA_hom/p.R).*(1/p.T-1/p.T_ref)); % Rate consatnts at reactor temperature(65 C) (m3/mol s)

        

        % ROTATION KINETICS
        p.G = p.h*p.rd^-1; % Normalized Gap distance 

        % INITIAL CONCENTRATION
        p.Cin = [p.C0_TG;....
                 0; ....
                 0; ....
                 0; ....
                 0; ....
                 0; ....
                 0; ....
                 0; ....
                 0; ....
                 0; ....
                 0; ....  
                 p.C0_MeOH;
                 ]; % Intial Bulk Concentrations (mol/m3)
       
        p.options=odeset('RelTol', 1e-8, 'AbsTol', 1e-8);  % required accuracy of the solution
        

    else 
    p=varargin{1};
 
    end
        
end
