function p = parameters_slurry(varargin)

    if nargin == 0
        % --- Shared Parameters ---
        p.R = 8.314; % J/mol.K
        p.g = 9.81;
        p.T_ref = 50 + 273.15;
        p.T = 55 + 273.15;
        p.ST = 22.5e-3; % Surface tesnison (N/m)

        % REACTOR and REACTION PARAMETERS
        p.h = 2e-3; % Gap distance in RS-SDR (m)
        p.rd = 65e-3; % spinner dia (m)
        p.G = p.h*p.rd^-1; % Normalized Gap distance 

        %% Physical Properties
        p.rho_TG = 918.8; % kg/m3
        p.rho_MeOH = 792; % kg/m3
        p.MW_TG = 794.326; % g/mol
        p.MW_MeOH = 32; % g/mol

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

        %% Catalyst Properties
        p.wt_cat = 0.3; % 30% wt of TG
        p.rho_cat = 2610; % kg/m3 (for Na2SiO3)
        p.dp = 100e-9; % particle diameter (m)
        p.eps_p = 0.5; % porocity of agglomerates
        p.tau_p = 3; % tortuosity of agglomerates
        p.L_sio = 42.8e3; % Bronsted sites (mol_SiO/m2_cat)
        p.Fv_s = p.wt_cat/(1-p.wt_cat)/p.rho_cat*(p.mdot_TG); % Volumetric flow rate of solid phase (m3_s/s)
        p.Deff = p.eps_p/p.tau_p*p.D_ME; % Diffusivites to the surface of the catalyst agglomerates (m2/s)

        %% Volume and Holdups
        p.V_r = 6.2e-5; % m3
        p.Fv_Tot = p.Fv_MeOH+p.Fv_TG+p.Fv_s; % total volumetric flow (m3_r/s)
        p.eps_TG = p.Fv_TG / p.Fv_Tot; % epsilon TG
        p.eps_MeOH = p.Fv_MeOH / p.Fv_Tot; % epsilon MeOH
        p.eps_s = p.Fv_s/p.Fv_Tot; % epsilon solid
        p.V_MeOH = p.eps_MeOH * p.V_r; % Volume of methanol (m3_MeOH)
        p.V_pores = p.eps_p*p.eps_s*p.V_r; % Volume of agglomarate pores (m3_p)
        p.V_solid = p.eps_s*p.V_r; % Voulme of the solid (m3_s)
        p.V_TG = p.eps_TG * p.V_r; % Volume of TG Phase (m3_TG)
        p.aLS = 6/p.dp*p.V_solid; % (m2_cat) aLS/V_r


        %% Kinetic Parameters (converted to SI)
        p.k_het_ref = [0.0842, 0.1273, 0.1032, 0.0232, 0.0571, 0.0083]/60/1000; % (m3/mol/s)/m2
        p.EA_het = [90.5e3, 87.3e3, 19.2e3, 2001.7e3, 28.3e3, 57.5e3]; % Activation energy (J/mol)
        p.k_het = p.k_het_ref .* exp((-p.EA_het/p.R).*(1/p.T-1/p.T_ref)); % Rate consatnts at reactor temperature(65 C) (m3/mol s)

        %% DISTRIBUTION COEFFICIENT
        p.m_DC = [150, 105, 54, 0.0062, 47, 0.0031]; % m3_MeOH/m3_TG (m_TG,m_DG,m_MG,m_G,m_FAME,m_MeOH)

        %% Initial Concentrations
        p.C0_TG = (p.rho_TG / p.MW_TG) * 1000;
        p.C0_MeOH = (p.rho_MeOH / p.MW_MeOH) * 1000;

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
                 0; ....
                 0; ....
                 0; ....
                 0; ....
                 0; ....
                 0; ....  
                 ]; % Intial Bulk Concentrations (mol/m3)

        p.options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

    else
        p = varargin{1};
    end
end
