% Code to Plot the variation of certain factor withe respect to the
% rotation speed and also to plot the concentration profile of one Rs_SDR

% Reset environment before simulation
clearvars;    % Clears all variables from workspace
close all;    % Closes all open figure windows
clc;          % Clears the command window

%% Accessing the inital variable in the structure
p = parameters_coated();

p.tspan = [0 60]; % Simulation Time(s)
rot_speed = 10:10:450; %rotation speeds(rad/s)

% Intializing Zero arrays for the analysis
E_DR =zeros(6,1); % To store the Dissipation energy for each rotational speed value
d_avg = zeros(6,1); % To store the average diameter of droplets for each rotational speed value
X = zeros(6,1); % To store the Conversion for each rotational speed value

for i = 1:length(rot_speed)

    %% General parameters calculations
    p.omega = rot_speed(i); % Rotational Speed (rad/s)
    p.Re_omg = p.omega * p.rd^2 /p.mu_k_TG; % Rotational Reynolds nummber
    p.Edr = (5.73*10^-12)*(p.G^-0.14)*(p.Re_omg^2.12); % Energy Disipation rate (w)
    E_DR(i) = p.Edr; % Storing Edr for each run
    p.spe_enrgy = p.Edr/(p.rho_TG*p.V_r); % Specific Energy Dissipation rate (W/kg)
    p.davg = 0.062*(p.ST/p.rho_TG)^(3/5)*p.spe_enrgy^(-2/5); % Average Droplet size for homogenous system for comparison (m)
    d_avg(i) = p.davg; % storing d_avg for each


    p.aLL = 6/p.davg*p.eps_MeOH; % Area of the droplet (m2) aLL/V_r
    % MASS TRANSFER COEFFICIENTS METHANOL PHASE
    p.Sh_MeOH = 2; % For Methanol Phase
    p.K_ME = (p.Sh_MeOH*p.D_ME)./p.davg; %  Mass transfer coefficient of components in ME (m3_i_ME/m2_int.s)

    % MASS TRANSFER COEFFICIENTS TG PHASE
    p.Sc = p.mu_k_TG ./ p.D_TG; % Schimdt number for all componenets in TG phase
    p.Sh_TG = 2+ (0.4*p.spe_enrgy*p.davg^4*p.mu_k_TG^-3)^(1/4)*p.Sc.^(1/3); % % Sherwood number for all componenets in TG phase
    p.K_TG = (p.Sh_TG.*p.D_TG)./p.davg; %  Mass transfer coefficient of components in TG (m3_i_TG/m2_int.s)

    % MASS TRANSFER COEFFICIENTS SOLID PHASE
    p.Sh_solid = 0.892*p.Re_omg^0.57*p.Sc.^(1/3) ; 
    
    % Overall MTC for CD interface
    p.K_OV_LL = 1./((1./p.K_TG) + (p.m_DC./p.K_ME)); % Overall Mass Transfer Coefficients (m3_i_TG/m2_int .s)

    % Overall MTC for DS interface
    p.K_OV_LS = (p.Sh_solid.*p.Deff)./p.L_coat; % Overall Mass Transfer Coefficients (m3_i_TG/m2_int .s)

    %solver
    [t, y, B] = ode15s(@(t, y) modelODE(t, y, p), p.tspan, p.Cin, p.options);
    fprintf('\nOil Phase: TG: %f, DG: %f, MG: %f, G: %f, FAME: %f, MeOH: %f\nMethanol Phase: TG: %f, DG: %f, MG: %f, G: %f, FAME: %f, MeOH: %f\nSolid Phase: TG: %f, DG: %f, MG: %f, G: %f, FAME: %f, MeOH: %f\n', ...
            y(end,1), y(end,2), y(end,3), y(end,4), y(end,5), y(end,6), y(end,7), y(end,8), y(end,9), y(end,10),y(end,11),y(end,12),y(end,13),y(end,14),y(end,15),y(end,16),y(end,17),y(end,18));
    p.C_end = [y(end,1);y(end,2);y(end,3);y(end,4);y(end,5);y(end,6);y(end,7);y(end,8);y(end,9);y(end,10);y(end,11);y(end,12);y(end,13);y(end,14);y(end,15);y(end,16);y(end,17);y(end,18);]; % End time concentration

    X(i) = ((p.F_TG) - ((p.C_end(7)*p.Fv_MeOH)+(p.C_end(1)*p.Fv_TG)))/(p.F_TG); %Storing conversion values

    P_FAME(i) = 3*X(i)*p.F_TG/p.V_r/60;%Storing productivity values
    

end

outputFolder = 'Plots';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end


%% Plots


%% conversion
%plotting the Effect of Rotation speed on conversion
figure;
plot(rot_speed,X, 'o--',color='red')
xlabel('Rotation speed (rad/s)')
ylabel('Conversion of TG (X)')
title('Effect of Rotation speed on conversion')
%saveas(gcf, fullfile(outputFolder, 'Conversion_vs_RotationSpeed.png'));

%% productivity
% plotting the Effect of Rotation speed on productivity
figure;
plot(rot_speed,P_FAME, '-',color='blue')
xlabel('Rotation speed (rad/s)')
ylabel('Productivity of FAMEs ($$\mathrm{mol_{FAME}}/\mathrm{m_R^3 \cdot min}$$)', ...
       'Interpreter', 'latex');
title('Effect of Rotation speed on productivity')
%saveas(gcf, fullfile(outputFolder, 'Productivity_vs_RotationSpeed.png'));

%% Energy dissipiation rate
% Plotting the Effect of Rotation speed on Energy Dissipation rate
figure;
plot(rot_speed,E_DR, 'o--',color='red')
xlabel('Rotation speed (rad/s)')
ylabel('Energy Dissipation rate (W)')
title('Effect of Rotation speed on Energy Dissipation rate ')
%saveas(gcf, fullfile(outputFolder, 'Edr_vs_RotationSpeed.png'));

% Plotting the Effect of Rotation speed on d_avg
figure;
plot(rot_speed,d_avg, 'o--',color='red')
xlabel('Rotation speed (rad/s)')
ylabel('Average droplet size (m)')
title('Effect of Rotation speed on davg ')
%saveas(gcf, fullfile(outputFolder, 'davg_vs_RotationSpeed.png'));

% Plotting the concentration profile in the Continuous (TG) phase only
legend_entries = {'TG_c', 'DG_c', 'MG_c', 'G_c', 'FAME_c', 'MeOH_c'};

figure;
plot(t, y(:,1:6)); % Columns 1 to 6 correspond to continuous phase
legend(legend_entries, 'Location', 'best');
xlabel('Time (s)');
ylabel('Concentration (mol/m^3_{TG})');
title('Concentration Profile in Continuous Phase (TG)');
grid on;
%saveas(gcf, fullfile(outputFolder, 'ConcentrationProfiles.png'));

% Store conversion for current run and catalyst thickness
%if exist('X_vs_cat_thick.mat', 'file')
%    load('X_vs_cat_thick.mat', 'cat_thickness', 'conversion_X');
%else
%    cat_thickness = [];
%    conversion_X = [];
%end

% Append current loading and conversion
%cat_thickness(end+1) = p.L_coat;
%conversion_X(end+1) = X(end); % or max(X) depending on what you want
%save('X_vs_cat_thick.mat', 'cat_thickness', 'conversion_X');

load('X_vs_cat_thick.mat', 'cat_thickness', 'conversion_X');

figure;
plot(cat_thickness, conversion_X, 'o--', 'LineWidth', 2);
xlabel('Catalyst thickness (m)');
ylabel('TG conversion');
title('Effect of Catalyst Thickness on Conversion');
grid on;

% Store conversion for current run and catalyst loading
%if exist('X_vs_T.mat', 'file')
    %load('X_vs_T.mat', 'Temp', 'conversion_X');
%else
    %Temp = [];
    %conversion_X = [];
%end

% Append current loading and conversion
%Temp(end+1) = p.T;
%conversion_X(end+1) = X(end); % or max(X) depending on what you want
%save('X_vs_T.mat', 'Temp', 'conversion_X');

% Plot Temp vs conversion
load('X_vs_T.mat', 'Temp', 'conversion_X');

figure;
plot(Temp, conversion_X, 'o--', 'LineWidth', 2);
xlabel('Temperature (K)');
ylabel('TG conversion');
title('Effect of Temperature on Conversion');
grid on;


%% Defining the ODE function
function [dydt, MT_LL, MT_LS] = modelODE(~, y, p)
    MT_LL_TG   = p.K_OV_LL(1)*p.aLL*(y(1)-p.m_DC(1)*y(7));
    MT_LL_DG   = p.K_OV_LL(2)*p.aLL*(y(2)-p.m_DC(2)*y(8));
    MT_LL_MG   = p.K_OV_LL(3)*p.aLL*(y(3)-p.m_DC(3)*y(9));
    MT_LL_G    = p.K_OV_LL(4)*p.aLL*(y(4)-p.m_DC(4)*y(10));
    MT_LL_FAME = p.K_OV_LL(5)*p.aLL*(y(5)-p.m_DC(5)*y(11));
    MT_LL_MeOH = p.K_OV_LL(6)*p.aLL*(y(6)-p.m_DC(6)*y(12));
    
    MT_LS_TG   = p.K_OV_LS(1)*p.aLS*(y(1)-y(13));
    MT_LS_DG   = p.K_OV_LS(2)*p.aLS*(y(2)-y(14));
    MT_LS_MG   = p.K_OV_LS(3)*p.aLS*(y(3)-y(15));
    MT_LS_G    = p.K_OV_LS(4)*p.aLS*(y(4)-y(16));
    MT_LS_FAME = p.K_OV_LS(5)*p.aLS*(y(5)-y(17));
    MT_LS_MeOH = p.K_OV_LS(6)*p.aLS*(y(6)-y(18));

    %R_1 = p.k_het(1)*y(13)*y(18)*p.aLS;   % TG + MeOH → FAME + DG
    %R_2 = p.k_het(2)*y(14)*y(18)*p.aLS;   % DG + MeOH → FAME + MG
    %R_3 = p.k_het(3)*y(15)*y(18)*p.aLS;   % MG + MeOH → FAME + G
    %R_4 = p.k_het(4)*y(14)*y(17)*p.aLS;   % FAME + DG → TG + MeOH
    %R_5 = p.k_het(5)*y(15)*y(17)*p.aLS;   % FAME + MG → DG + MeOH
    %R_6 = p.k_het(6)*y(16)*y(17)*p.aLS;   % FAME + G → MG + MeOH

    %TG phase balance
    dNTG_o_dt= (p.Fv_TG*p.Cin(1)-p.Fv_TG*y(1)-MT_LL_TG-MT_LS_TG)/p.V_MeOH; 

    dNDG_o_dt= (                -p.Fv_TG*y(2)-MT_LL_DG-MT_LS_DG)/p.V_MeOH;  

    dNMG_o_dt= (                -p.Fv_TG*y(3)-MT_LL_MG-MT_LS_MG)/p.V_MeOH;

    dNG_o_dt=  (                -p.Fv_TG*y(4)-MT_LL_G-MT_LS_G)/p.V_MeOH;

    dNF_o_dt=  (                -p.Fv_TG*y(5)-MT_LL_FAME-MT_LS_FAME)/p.V_MeOH;

    dNME_o_dt= (                -p.Fv_TG*y(6)-MT_LL_MeOH-MT_LS_MeOH)/p.V_MeOH;

    % Methanol phase

    dNTG_m_dt= (                   -p.Fv_MeOH*y(7) +MT_LL_TG)/p.V_TG;

    dNDG_m_dt= (                   -p.Fv_MeOH*y(8) +MT_LL_DG)/p.V_TG;
   
    dNMG_m_dt= (                   -p.Fv_MeOH*y(9) +MT_LL_MG)/p.V_TG; 
    
    dNG_m_dt=  (                   -p.Fv_MeOH*y(10)+MT_LL_G)/p.V_TG;
    
    dNF_m_dt=  (                   -p.Fv_MeOH*y(11)+MT_LL_FAME)/p.V_TG;
    
    dNME_m_dt= (p.Fv_MeOH*p.Cin(12)-p.Fv_MeOH*y(12)+MT_LL_MeOH)/p.V_TG;

    % Solid phase

    dNTG_s_dt=( +MT_LS_TG - p.k_het(1)*y(13)*y(18)*p.aLS + p.k_het(4)*y(14)*y(17)*p.aLS)/p.V_solid;

    dNDG_s_dt= (+MT_LS_DG + p.k_het(1)*y(13)*y(18)*p.aLS - p.k_het(2)*y(14)*y(18)*p.aLS - p.k_het(4)*y(14)*y(17)*p.aLS + p.k_het(5)*y(15)*y(17)*p.aLS)/p.V_solid;
    
    dNMG_s_dt= ( +MT_LS_MG + p.k_het(2)*y(14)*y(18)*p.aLS - p.k_het(3)*y(15)*y(18)*p.aLS - p.k_het(5)*y(15)*y(17)*p.aLS + p.k_het(6)*y(16)*y(17)*p.aLS)/p.V_solid; 
    
    dNG_s_dt= ( +MT_LS_G + p.k_het(3)*y(15)*y(18)*p.aLS - p.k_het(6)*y(16)*y(17)*p.aLS)/p.V_solid;
    
    dNF_s_dt= (  +MT_LS_FAME + p.k_het(1)*y(13)*y(18)*p.aLS + p.k_het(2)*y(14)*y(18)*p.aLS + p.k_het(3)*y(15)*y(18)*p.aLS - p.k_het(4)*y(14)*y(17)*p.aLS -p.k_het(5)*y(15)*y(17)*p.aLS - p.k_het(6)*y(16)*y(17)*p.aLS)/p.V_solid;
    
    dNME_s_dt=( +MT_LS_MeOH - p.k_het(1)*y(13)*y(18)*p.aLS - p.k_het(2)*y(14)*y(18)*p.aLS - p.k_het(3)*y(15)*y(18)*p.aLS + p.k_het(4)*y(14)*y(17)*p.aLS + p.k_het(5)*y(15)*y(17)*p.aLS + p.k_het(6)*y(16)*y(17)*p.aLS)/p.V_solid;


    dydt = [dNTG_o_dt; dNDG_o_dt; dNMG_o_dt; dNG_o_dt; dNF_o_dt; dNME_o_dt; dNTG_m_dt; dNDG_m_dt; dNMG_m_dt; dNG_m_dt; dNF_m_dt; dNME_m_dt; dNTG_s_dt; dNDG_s_dt; dNMG_s_dt; dNG_s_dt; dNF_s_dt; dNME_s_dt;];
    MT_LL = [MT_LL_TG; MT_LL_DG; MT_LL_MG; MT_LL_G; MT_LL_FAME; MT_LL_MeOH];
    MT_LS = [MT_LS_TG; MT_LS_DG; MT_LS_MG; MT_LS_G; MT_LS_FAME; MT_LS_MeOH];
end


