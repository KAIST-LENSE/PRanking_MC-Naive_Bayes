%% Modelling of Biomass Growth inside Vertical Bubble Column Cultivators
% Primary Reference: Tebbani S. CO2 Biofixation by Microalgae: Modeling, 
%                    Estimation and Control. ISTE, 2014
function dX = Cultivation_Model(~,Input)
%% Define Model Inputs
% Input(1) = X                  % Concentration of 10^9 cells per Liter
% Input(2) = TIC                % Total Inorganic Carbon, mmol/L
% Variable(3) = Mu              % Specific Biomass Growth Rate, 1/hr
% Variable(4) = CO2             % CO2 concentration in media, mol/L
% Variable(5) = E               % Average Irradiance per Column
% Variable(6) = I_out           % Outgoing light as a fx of incident light
% Import Global Parameters
global K_CL K_1 K_2 K_E kLa F_CUL CO2e H_ion Mu_max V_col Yc A_r I_in C1 C2

%% Cultivation Model Equations
% CO2 concentration in media as a function of pH and dissociation
CO2 = Input(2)/(1 + (K_1/H_ion) + ((K_1*K_2)/H_ion^2));
% Outgoing Light Intensity
I_out = C1*I_in*Input(1)^C2;
% Average Irradiance per Column
E = ((I_in - I_out)*A_r)/(V_col*Input(1));
% Specific Growth Rate via Monod Kinetics, limited by [TIC]
Mu = Mu_max*(E/(K_E + E))*(Input(2)/((K_CL*Input(1)) + Input(2)));
% Change in Biomass Concentration due to Growth
dX(1) = Input(1)*Mu - (F_CUL/V_col)*Input(1);
% Change in [TIC] due to Biomass fixation, under constant CO2 feed
dX(2) = -(Input(1)*Mu/Yc) + kLa*(CO2e - CO2);

%% Export Differential Equations
% Equations for Specific Growth Rate and CO2 concentration are Implicit
dX = [dX(1); dX(2)];
end