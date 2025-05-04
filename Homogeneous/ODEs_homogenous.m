% Code to Plot the variation of certain factor withe respect to the
% rotation speed and also to plot the concentration profile of one Rs_SDR

% Reset environment before simulation
clearvars;    % Clears all variables from workspace
close all;    % Closes all open figure windows
clc;          % Clears the command window

%% Accessing the inital variable in the structure
p = parameters_homogeneous();

p.tspan = [0 120]; % Simulation Time(s)
rot_speed = 10:10:450; %rotation speeds(rad/s)

% Intializing Zero arrays for the analysis
E_DR =zeros(6,1); % To store the Dissipation energy for each rotational speed value
d_avg = zeros(6,1); % To store the average diameter of droplets for each rotational speed value
X = zeros(6,1); % To store the Conversion for each rotational speed value


for i = 1:length(rot_speed)
    p.omega = rot_speed(i); % Rotational Speed (rad/s)
    p.Re_omg = p.omega * p.rd^2 /p.mu_k_TG; % Rotational Reynolds nummber
    p.Edr = (5.73*10^-12)*(p.G^-0.14)*(p.Re_omg^2.12); % Energy Disipation rate (w)
    E_DR(i) = p.Edr; % Storing Edr for each run
    p.spe_enrgy = p.Edr/(p.rho_TG*p.V_r); % Specific Energy Dissipation rate (W/kg)
    p.davg = 0.062*(p.ST/p.rho_TG)^(3/5)*p.spe_enrgy^(-2/5); % Average Droplet size (m)
    d_avg(i) = p.davg; % storing d_avg for each
    p.Sa = 6/p.davg; % Specific surface area of the droplet (m2/m3_MeOH)

    % MASS TRANSFER COEFFICIENTS METHANOL PHASE
    p.Sh_MeOH = 2; % For Methanol Phase
    p.K_ME = (p.Sh_MeOH*p.D_ME)./p.davg; %  Mass transfer coefficient of components in ME (m3_i_ME/m2_int.s)

    % MASS TRANSFER COEFFICIENTS TG PHASE
    p.Sc = p.mu_k_TG ./ p.D_TG; % Schimdt number for all componenets in TG phase
    p.Sh_TG = 2+ (0.4*p.spe_enrgy*p.davg^4*p.mu_k_TG^-3)^(1/4)*p.Sc.^(1/3); % % Sherwood number for all componenets in TG phase
    p.K_TG = (p.Sh_TG.*p.D_TG)./p.davg; %  Mass transfer coefficient of components in TG (m3_i_TG/m2_int.s)

    % Overall MTC
    p.K_OV = 1./((1./p.K_TG) + (p.m_DC./p.K_ME)); % Overall Mass Transfer Coefficients (m3_i_TG/m2_int .s)

    %solver
    [t, y, B] = ode15s(@(t, y) modelODE(t, y, p), p.tspan, p.Cin, p.options);
    fprintf('\nOil Phase: TG: %f, DG: %f, MG: %f, G: %f, FAME: %f, MeOH: %f\nMethanol Phase: TG: %f, DG: %f, MG: %f, G: %f, FAME: %f, MeOH: %f\n', ...
            y(end,1), y(end,2), y(end,3), y(end,4), y(end,5), y(end,6), y(end,7), y(end,8), y(end,9), y(end,10),y(end,11),y(end,12));
    p.C_end = [y(end,1);y(end,2);y(end,3);y(end,4);y(end,5);y(end,6);y(end,7);y(end,8);y(end,9);y(end,10);y(end,11);y(end,12);]; % End time concentration

    X(i) = ((p.F_TG) - ((p.C_end(7)*p.Fv_MeOH)+(p.C_end(1)*p.Fv_TG)))/(p.F_TG); %Storing conversion values

    P_FAME(i) = 3*X(i)*p.F_TG/p.V_r/60;%Storing productivity values

end

outputFolder = 'Plots';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end


%% Plots


%% conversion
% plotting the Effect of Rotation speed on conversion
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

% Store conversion for current run and catalyst loading
%if exist('X_vs_T.mat', 'file')
%    load('X_vs_T.mat', 'Temp', 'conversion_X');
%else
%    Temp = [];
%    conversion_X = [];
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
function [dydt, MT] = modelODE(~, y, p)
    MT_TG = p.K_OV(1)*p.Sa*(y(1)-p.m_DC(1)*y(7));
    MT_DG =p.K_OV(2)*p.Sa*(y(2)-p.m_DC(2)*y(8));
    MT_MG = p.K_OV(3)*p.Sa*(y(3)-p.m_DC(3)*y(9));
    MT_G = p.K_OV(4)*p.Sa*(y(4)-p.m_DC(4)*y(10));
    MT_FAME = p.K_OV(5)*p.Sa*(y(5)-p.m_DC(5)*y(11));
    MT_MeOH = p.K_OV(6)*p.Sa*(y(6)-p.m_DC(6)*y(12));
   

    %TG phase balance
    dNTG_o_dt=((p.Fv_TG*p.Cin(1)-p.Fv_TG*y(1)-MT_TG))/p.V_TG; 

    dNDG_o_dt=((-MT_DG-p.Fv_TG*y(2)))/p.V_TG;  

    dNMG_o_dt=((-MT_MG-p.Fv_TG*y(3)))/p.V_TG;  

    dNG_o_dt=((-MT_G-p.Fv_TG*y(4)))/p.V_TG; 

    dNF_o_dt=((-MT_FAME-p.Fv_TG*y(5)))/p.V_TG;

    dNME_o_dt=((-MT_MeOH - p.Fv_TG*y(6)))/p.V_TG;

    % Methanol phase

    dNTG_m_dt=((MT_TG-p.k_hom(1)*y(7)*y(12)*p.V_MeOH+p.k_hom(4)*y(11)*y(8)*p.V_MeOH-p.Fv_MeOH*p.eps_MeOH*y(7)))/p.V_MeOH; 

    dNDG_m_dt=((MT_DG + p.k_hom(1)*y(7)*y(12)*p.V_MeOH-p.k_hom(2)*y(8)*y(12)*p.V_MeOH+p.k_hom(5)*y(11)*y(9)*p.V_MeOH-p.k_hom(4)*y(8)*y(11)*p.V_MeOH - p.Fv_MeOH*p.eps_MeOH*y(8) ))/p.V_MeOH;
    
    dNMG_m_dt=((MT_MG+p.k_hom(2)*y(8)*y(12)*p.V_MeOH-p.k_hom(3)*y(9)*y(12)*p.V_MeOH-p.k_hom(5)*y(9)*y(11)*p.V_MeOH+p.k_hom(6)*y(11)*y(10)*p.V_MeOH - p.Fv_MeOH*p.eps_MeOH*y(9) ))/p.V_MeOH; 
    
    dNG_m_dt=((MT_G+p.k_hom(3)*y(9)*y(12)*p.V_MeOH-p.k_hom(6)*y(11)*y(10)*p.V_MeOH - p.Fv_MeOH*p.eps_MeOH*y(10)))/p.V_MeOH;
    
    dNF_m_dt=((MT_FAME +p.k_hom(1)*y(7)*y(12)*p.V_MeOH+p.k_hom(2)*y(8)*y(12)*p.V_MeOH+p.k_hom(3)*y(9)*y(12)*p.V_MeOH-p.k_hom(4)*y(8)*y(11)*p.V_MeOH-p.k_hom(5)*y(9)*y(11)*p.V_MeOH-p.k_hom(6)*y(11)*y(10)*p.V_MeOH - p.Fv_MeOH*p.eps_MeOH*y(11)))/p.V_MeOH;
    
    dNME_m_dt=((p.Fv_MeOH*p.Cin(12)-p.Fv_MeOH*p.eps_MeOH*y(12)-MT_MeOH-p.k_hom(1)*y(7)*y(12)*p.V_MeOH-p.k_hom(2)*y(8)*y(12)*p.V_MeOH-p.k_hom(3)*y(9)*y(12)*p.V_MeOH+p.k_hom(4)*y(8)*y(11)*p.V_MeOH+p.k_hom(5)*y(9)*y(11)*p.V_MeOH+p.k_hom(6)*y(11)*y(10)*p.V_MeOH))/p.V_MeOH;



    dydt = [dNTG_o_dt; dNDG_o_dt; dNMG_o_dt; dNG_o_dt; dNF_o_dt; dNME_o_dt; dNTG_m_dt; dNDG_m_dt; dNMG_m_dt; dNG_m_dt; dNF_m_dt; dNME_m_dt;];
    MT = [MT_TG; MT_DG; MT_MG; MT_G; MT_FAME; MT_MeOH];
end


