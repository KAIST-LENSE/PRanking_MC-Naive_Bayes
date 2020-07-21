%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CCU Model for Microalgal Biomass-based Biocrude Production Process %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EvalMetrics = CCUS_Biocrude(kp)             
%% Declare Global Variables
global K_CL K_1 K_2 K_E kLa F_CUL CO2e H_ion Mu_max V_col Yc A_r I_in C1 C2
%% KEY PARAMETERS WITH UNCERTAINTY
% [1] Cultivation Process
avg_GPC = kp(1);       % Average dry weight, 0.0224 grams/10^9 cells
Mu_max = kp(2);        % Maximum Specific Growth Rate, 1.07 hr^-1
kLa = kp(3);           % Volumetric mass transfer coeff, 1.4 1/hr
y_Lipid = kp(4);       % PSEUDO Cell Lipid 30% to match Oil Yields   
P_CO2 = kp(5);         % Partial Pressure of CO2 from FG feed, 0.1358 atm
Daylight = kp(6);      % Daily fraction available for photosynthesis

% [2] Harvest+Dewater Process
rec_CENTR = kp(7);     % Biomass recovery during centrifuge harvest, 95%

% [3] Blowdown + Media Recycle
% [4] Wet Biomass HTL Treatment
% [5] Crude Stream Refining Process
% [6] Oil Stream Hydrotreating Process
% [TEA] Techno-Economic Analysis Parameters
C_PBR = kp(8);         % 65153 $/Acre

% [LCA] CO2 Life Cycle Assessment Parameters
GWI_ELEC_GRID = kp(9); % 52.78 kg CO2-eq/GJ Electricity
GWI_NSOURCE = kp(10);  % 3.01 kg CO2-eq/kg N-Nutrients



%% NON VARYING FIXED-PARAMETERS & DESIGN SPECS
% [1] Cultivation Process
pH = 9.5;              % Alkaline pH control setpoint for optimal growth
V_col = 6.25;          % Design Spec Vol. of Cult Bag (KCRC=6.25L)
K_CL = 3.8;            % Half-TIC sat. constant, mmol/10^9 cells
K_1 = 10^-6.35;        % Dissociation constant of CO2, mol/L
K_2 = 10^-10.3;        % Dissociation constant of Bicarbonate, mol/L
K_E = 0.08;            % Light half sat. constant, mu*E/s
I_in = 90;             % Average incident light intensity, mu*E/m^2*s
A_r = 0.31;            % Illuminated area per unit PBR, m^2
C1 = 0.493;            % Reactor geometry constant
C2 = -0.925;           % Light wavelength constant
Yc = 1211;             % Conversion Yield of Carbon, in 10^9 cells
F_CUL = 0.06;          % Medium cycle flowrate (REF=0.12L/hr,KCRC=0.06L/hr)
VVM = 0.05;            % Flowrate of flue gas in VVM (from Korea Univ.)
RTP_FG = 0.0264385;    % RT/P of Flue Gas at 45 deg C
H = 29.41;             % Henry's Constant, atmL/mol
X_inoc = 3.12;         % Average innoculum concentration in 10^9 cells/L
TIC_0 = 0.0046;        % Initial TIC concentration in media, moles/L
HAR_T = 96;            % Semi-continuous harvest time in hours, 96 hrs 
STOCK_T = 0.5;         % Time to transfer stock soln. to GM Tank, hrs
GM_T = 1.0;            % Time to transfer growth media to Cult., hrs
REC_T = 1.0;           % Time to transfer cult. broth back to GM Tank, hrs
V_Cult_T = 1000000000; % Total cultivation volume, liters (1E6 m^3)
y_Prote = 0.22;        % PSEUDO Cell Protein % to match HTL Oil Yields
y_Ash = 0.03;          % PSEUDO Cell Ash % to match HTL Oil Yields
V_Cult_Unit = 36.3;    % m^3 per unit Cultivation Module
Land_Cult_Unit = 1250; % m^2 per unit Cultivation Module
FG_MM = [44.01, 18.02, 34.1, 100.09, 46.01, 32, 28.01]; 
FG_Moles = [P_CO2, 0.0818, 0.0005, 0.0154, 0.0025, 0.0354, 0.7286];

% [2] Harvest+Dewater Process
conc_CENTR = 0.25;     % Exit conc. of biomass @ Centrifuge, 25 dwt%
SEP_T = 10;            % Operation time/batch for culture separation, hrs

% [3] Blowdown + Media Recycle
% [4] Wet Biomass HTL Treatment
HTL_Aspen = load('HTL_MassBalance_Data');% Load HTL Reactor MB Data 
HTL_Sto = HTL_Aspen.Stoic_Coeff;         % Stoic coefficients for HTL rxns
HTL_Key = HTL_Aspen.Rxn_Key_Comp;        % Key reactant matrix
HTL_Exit = HTL_Aspen.Exit_Frac;          % Phase sep. partition coefficient
X_C = 0.85;            % Conversion of cell carbohydrates
X_L = 0.6;             % Fractional conversion of biomass lipid in HTL, 60%
X_P1 = 0.65;           % Conversion of cell proteins (Type A)
X_P2 = 0.65;           % Conversion of cell proteins (Type B)
X_P3 = 0.65;           % Conversion of cell proteins (Type C)
X_P4 = 0.65;           % Conversion of cell proteins (Type D)
X_NH3 = 0.70;          % Secondary conversion of NH3
X_AA = 0.68;           % Secondary reduction of Acetic Acids
X_LA = 0.85;           % Secondary reduction of Lactic Acids
X_GL = 0.85;           % Secondary reduction + breakdown of glycerols
X_FA = 0.69;           % Secondary reduction + breakdown of Fatty Acids
delP_Slry = 138;                         % Pressure increase (bar) for HTL
SG_MA = 1.07;                            % Specific Gravity 1/kg, [18]
BD_LOSS = [0.01, 0.01, 0, 0, 0];         % Blowdown Loss Fractions

% [5] Crude Stream Refining Process
REF_Aspen = load('REF_MassBalance_Data');% Load Oil-Phase Refining MB Data
REF_Frac = REF_Aspen.Separation_Frac;

% [6] Oil Stream Hydrotreating Process
UPG_Aspen = load('UPG_MassBalance_Data');% Load Crude Upgrading MB Data
UPG_Sto = UPG_Aspen.Stoic_Coeff;         % Stoich coefficients for upg rxns
UPG_Key = UPG_Aspen.Rxn_Key_Comp;        % Key reactant matrix
UPG_FC = UPG_Aspen.Frac_Conv;            % Fractional Conversion matrix
UPG_Exit = UPG_Aspen.Exit_Frac;          % Fractionation separation coeffs

% [TEA] Techno-Economic Analysis Parameters
%%% Chemical Engineering Plant Cost Indices (CEPCI)
CEPCI = 619.2;                        % Based on CEPCI for 2019
CEPCI_PBR = 619.2;                    % 2019 CEPCI based on [21] publish yr
CEPCI_PumpsBlowers = 541.3;           % 2016 CEPCI based on [15] publish yr
CEPCI_Separator = 541.3;              % 2016 CEPCI based on [15] publish yr
CEPCI_Upstr_Tanks = 603.1;            % 2018 CEPCI from [22]
CEPCI_DS = 567.5;                     % 2017 CEPCI for Downstream (APEAv10)
%%% Equipment Efficiencies
NU_Pump = 0.73;                       % Efficiency for Water & Solutions
NU_Motor = 0.90;                      % Motor efficiency for 50-250kW Range
%%% Nth-Plant Assumption Exponential Constants
NTH_PMP = 0.6;                        % From [23]
NTH_HEX = 1.2;                        % From [23], Shell and Tube
NTH_CSTR = 0.8;                       % From [23], Jacketed CSTR
NTH_DIST = 0.85;                      % From [23], Vertical Pressurized
NTH_COMP = 0.8;                       % From [23], Blower Compressors
%%% Capital Cost Lang Factors
OSBL_OS = 0.4;                        % Offsite Operations, from [23]
OSBL_DE = 0.25;                       % Design and Engineering, from [23]
OSBL_CN = 0.1;                        % Contingency, from [23]
LIFET = 30;                           % Plant Lifetime, years
DISCO = 0.07;                         % Annual Discount Rate % of TCI
%%% [CAPEX] Quoted Equipment Cost Data (Pre-Sized Installed)
C_FG_Blower = 5803.2;                 % 7.5kW 3-lobe Blower, [15]
C_Stock_Pump = 1826.4;                % 5.5kW Centrifugal Pump, [15]
C_GM_Pump = 20234.7;                  % 10kW Submersible Pump, [15]
C_SEP = 62591.2;                      % 7.5kW Centrifugal Separator
C_Rec_Pump = 1826.4;                  % 5.5kW Centrifugal Pump, [15]
%%% [CAPEX] Equipment Cost Data (Nth-Plant Sizing, Installed)
REFCOST_SS_T = 54159;                 % Cost per Tank, Carbon Steel
REFCOST_GM_T = 478113;                % Cost per Tank, Carbon Steel
REFCOST_P100 = 889400;                % Ref. Bare Equipment Cost
REFCOST_H100 = 2479500;               % Ref. Bare Equipment Cost
REFCOST_E100 = 52000000;              % Based on 25333 kgDCW/hr from [24] 
REFCOST_H300 = 97200;                 % Ref. Bare Equipment Cost
REFCOST_C100 = 555500;                % Ref. Bare Equipment Cost
REFCOST_C200 = 283700;                % Ref. Bare Equipment Cost
REFCOST_C300 = 307900;                % Ref. Bare Equipment Cost
REFCOST_P300 = 209500;                % Ref. Bare Equipment Cost
REFCOST_P400 = 1676400;               % Ref. Bare Equipment Cost
REFCOST_F100 = 302200;                % Ref. Bare Equipment Cost
%%% [CAPEX] Equipment Size Data (Nth-Plant Sizing, Installed)
CAP_SS_T = 500;                       % m^3 Storage Capacity for SS Tank
CAP_GM_T = 10000;                     % m^3 Storage Capacity for GM Tank
%%% [OPERATING COST] Plant Costs
OP_PE_LAND = 340;        % $/Acre from [21]. Land maintainence costs
OP_PE_CIP = 750;         % $/Acre from [21]. Clean-in-place costs 
OP_PE_PBR = 5417;        % $/Acre from [21]. PBR plastic costs 6yr lifetime
%%% [OPERATING COST] Raw Material Costs
OP_RM_NSOURCE = 0.909;   % $/kg of N-Nutrients (Assumed Ammonia, [21])
OP_RM_PSOURCE = 0.742;   % $/kg of P-Nutrients (Assumed DAP, [21])
OP_RM_SSOURCE = 0.07;    % $/kg, Industry Quote Metal Sulfate fertilizers
OP_RM_FG = 0;                         % $/kg of Flue Gas
OP_RM_H2O = 0;                        % [21]. Depends on site   
OP_RM_SOLV = 0.495;                   % $/kg from Alibaba
OP_RM_H2 = 1.5;                       % $/kg from [25]
%%% [OPERATING COST] Utility Costs
UT_ELEC = 27.21;                      % $/GJ
UT_MP_STEAM = 5.303;                  % $/GJ based on 50 psig
UT_HP_STEAM = 5.194;                  % $/GJ based on 200 psig
UT_CW = -5.34;                        % $/GJ from [21]

% [LCA] CO2 Life Cycle Assessment Parameters
%%% [BOUNDARY] Binary Parameters to Set LCA Boundary for Plant
B_DIR = 1;    % Include direct plant emissions? Default=1
B_ENR = 1;    % '' indirect effects from energy consumption? Default=1
B_MAT = 1;    % '' indirect effects from RawMat. consumption? Default=1
B_CaS = 0;    % '' indirect effects from plant construction & salvage?
B_PC = 0;     % '' indirect effects from product consumption?
%%% [UTILITY] Global Warming Potential Parameters for Utility Consumption
GWI_STEAM_HP = 203.61;                % kg CO2-eq/GJ HP Steam
GWI_STEAM_MP = 222.25;                % kg CO2-eq/GJ MP Steam
GWI_COOLW = 0;                        % kg CO2-eq/GJ Cooling Water 
%%% [MATERIAL] Global Warming Potential Parameters for Material Consumption
GWI_CO2_FG = -1;                      % kg CO2-eq/kg CO2 in FG
GWI_PSOURCE = 1.46;                   % kg CO2-eq/kg P-Nutrients
GWI_SSOURCE = 0.3;                    % kg CO2-eq/kg S-Nutrients 
GWI_WATER = 0;                        % kg CO2-eq/kg Water Consum.
GWI_SOLV = 2.79;                      % kg CO2-eq/kg Decane Solv.
GWI_H2 = 1.63;                        % kg CO2-eq/kg H2 (liq)



%% DERIVED PARAMETERS AND DESIGN SPECS
% [1] Cultivation Process
CO2e = P_CO2/H;         % CO2 conc. (moles/L) that is in eq. with gas phase
H_ion = 10^-pH;         % Concentration of Hydrogen Ions in medium
N_col = V_Cult_T/V_col; % Number of Plastic Columns in Cultivation Plant
y_Carb = 1-y_Lipid-y_Prote-y_Ash;   % PSEUDO Cell Carbohydrate Content
% [2] Harvest+Dewater Process
% [3] Blowdown Recycle Treatment
% [4] Wet Biomass HTL Treatment
% [5] Crude Stream Refining Process
% [6] Oil Stream Hydrotreating Process
% [TEA] Techno-Economic Analysis Parameters
OP_Eff = HAR_T/(HAR_T+STOCK_T+SEP_T+max([GM_T,REC_T])); % Operating eff, %
Num_Har = floor((OP_Eff*8760)/HAR_T);                   % # of Harvests/yr
% [LCA] CO2 Life Cycle Assessment Parameters



%% PROCESS MASS BALANCE MODEL
% [1] Cultivation Process
% (a) Solve Cultivation Mass Balance Model
NUM_INT = HAR_T*2;                          % Number of Timespan Intervals
timespan = linspace(0, HAR_T, NUM_INT);
cult_initial = [X_inoc; TIC_0];             % Initial [TIC] is 0.0046 mol/L
[~,cult_out] = ode45(@Cultivation_Model, timespan, cult_initial);
% (b) Calculate Biomass Production
%%%Calculate Cell Concentration of the PBRs at each timestep
CELL_CONC = zeros(NUM_INT,1);
for i = 1:NUM_INT
    if i == 1 
        CELL_CONC(i) = X_inoc + cult_out(i,1);
    else
        CELL_CONC(i) = CELL_CONC(i-1) + cult_out(i,1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CULTIVATION ANALYTICS
% Plot Change in Cell Concentration over Time
%plot(timespan,CELL_CONC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Har_Conc = CELL_CONC(NUM_INT)*avg_GPC*(1/V_col);     % g/L conc. @ harvest
Inn_Conc = X_inoc*avg_GPC*(1/V_col);                 % g/L conc. @ innocul.
prod_bm = (Har_Conc - Inn_Conc)*V_col*N_col/1000;    % kg biomass @ harvest
% (c) Calculate Fraction of Cultivation that is Harvested
Har_Frac = 1-(Inn_Conc/Har_Conc);
% (d) Complete Cultivation Mass Balance
%%%Initial Nutrient loading of [N, P, S], (kg/L)
nutri_conc = [0.250 0.250 0.0904];                
%%%Specific (kg/kg) Biomass Consumption of [Water, CO2, N, P, S]
cult_stoic = [0.829, 2.037, 0.238, 0.009, 0.008];   
%%%Cultivation Raw Material Consumption
cons_cult = cult_stoic.*prod_bm;                    
% (e) Exit Flow Composition: [Water, Biomass, N, P, S], kg/hr 
exit_cult = [];
exit_cult(1) = (V_Cult_T - cons_cult(1))*Har_Frac./HAR_T;
exit_cult(2) = prod_bm*Har_Frac/HAR_T;
exit_cult(3:5) = ((V_Cult_T*nutri_conc)-cons_cult(3:5)).*(Har_Frac/HAR_T);

% [2] Harvest+Dewater Process 
% (a) Centrifugation
bm_CENTR = rec_CENTR*exit_cult(2);             % Recovered biomass, kg/hr
water_CENTR = (bm_CENTR/conc_CENTR)-bm_CENTR;  % Exit water flowrate, kg/hr
exit_CENTR = [water_CENTR, bm_CENTR, water_CENTR*nutri_conc]; 
% (b) Recycle Culture to Blowdown
recycle_CENTR = exit_cult - exit_CENTR;        % kg/hr

% [3] Blowdown + Media Recycle
% (a) Exit Streams
blowdown_loss = recycle_CENTR.*BD_LOSS;                   % kg/hr
cult_recycle = recycle_CENTR - blowdown_loss;             % kg/hr
wet_bm = exit_CENTR(1:2);                                 % kg/hr 
% (b) Make-Up Stream, kg
MU_Nutri = (nutri_conc.*V_Cult_T*Har_Frac) - (cult_recycle(3:5).*HAR_T);
MU_Water = (V_Cult_T*Har_Frac) - (cult_recycle(1).*HAR_T);

% [4] Wet Biomass HTL Treatment
% (a) Biomass Molar Breakdown [Water, Carbo, Lipids, Ash, Proteins], kg/hr
comp_biom = [exit_CENTR(1), exit_CENTR(2)*y_Carb, exit_CENTR(2)*y_Lipid,... 
             exit_CENTR(2)*y_Ash, exit_CENTR(2)*y_Prote*0.25,...
             exit_CENTR(2)*y_Prote*0.25, exit_CENTR(2)*y_Prote*0.25,...
             exit_CENTR(2)*y_Prote*0.25];
molarmass_biom = [18.02, 180.1, 805.4, 56.07, 131.1, 117.1, 89.09, 105];
%%%Molar Flowrate into HTL [Water, Carbo, Lipid, Ash, Proteins], kmole/hr
stream_biom = comp_biom./molarmass_biom;
% (b) Hydrothermal Liquefaction Reactor, kmole/hr 
%%%Inlet molar flow of Biomass Components into HTL reactor
HTL_In = zeros(1,34);
HTL_In(1:8) = stream_biom;     % Non-biomass components have 0 initial feed
%%%Define Fractional Conversions for Biomass Components during HTL
HTL_FC = [X_C; X_L; X_P1; X_P2; X_P3; X_P4; X_NH3; X_AA; X_LA; X_GL; X_FA];
%%%Extraction Solvent Addition Stream, kmole/hr
ExSolv = 0.481*wet_bm(2)/142.29;
%%%Calculate HTL Outlet Flow via Mass Balance Calculation, kmole/hr
Num_HTL_Rxn = size(HTL_Sto,1);
HTL_Out = [];
for i = 1:Num_HTL_Rxn
    if i == 1
        HTL_Out = HTL_In + sum((HTL_FC(i)*HTL_Key(i,:)).*HTL_In).*HTL_Sto(i,:);
    else
        HTL_Out = HTL_Out + sum((HTL_FC(i)*HTL_Key(i,:)).*HTL_Out).*HTL_Sto(i,:);
    end
end
%%%Addition of Oil Phase Extraction Solvent, kmole/hr
HTL_Out = [HTL_Out, ExSolv];
%%%HTL Process Exit Streams, kmole/hr
exit_HTLGas = HTL_Out.*HTL_Exit(1,:);
exit_HTLOil = HTL_Out.*HTL_Exit(2,:);
exit_HTLAq = HTL_Out.*HTL_Exit(3,:);
exit_HTLSol = HTL_Out.*HTL_Exit(4,:);

% [5] Crude Stream Refining Process
%%%Refining Process Separation Mass Balances, kmole/hr
C100_Top = exit_HTLOil.*REF_Frac(1,:);
C100_Bot = exit_HTLOil - C100_Top;
C200_Top = C100_Top.*REF_Frac(2,:);
C200_Bot = C100_Top - C200_Top;
C300_Top = C200_Bot.*REF_Frac(3,:);
%%%Consolidate Biocrude Stream, kmole/hr
exit_Biocrude = C300_Top + C100_Bot;

% [6] Oil Stream Hydrotreating Process
%%%Hydrogen Addition Stream, kmole/hr
UPG_In = exit_Biocrude;
UPG_In(36:40) = 0;
UPG_In(14) = 2.073*(exit_Biocrude(3) + exit_Biocrude(9) + ...
             sum(exit_Biocrude(17:18)) + exit_Biocrude(32));
%%%Calculate Upgraded Outlet Flow via Mass Balance Calculation, kmole/hr
Num_UPG_Rxn = size(UPG_Sto,1);
UPG_Out = [];
for i = 1:Num_UPG_Rxn
    if i == 1
        UPG_Out = UPG_In + sum((UPG_FC(i)*UPG_Key(i,:)).*UPG_In).*UPG_Sto(i,:);
    else
        UPG_Out = UPG_Out + sum((UPG_FC(i)*UPG_Key(i,:)).*UPG_Out).*UPG_Sto(i,:);
    end
end
%%%Fractionation Exit Streams, kmole/hr
exit_H2Rec = UPG_Out.*UPG_Exit(1,:);
exit_PROD = UPG_Out.*UPG_Exit(2,:);
MMass = load('Molar_Mass');
molar_mass = MMass.molmass;

% ORGANIZE PRODUCT STREAMS
PROD_BIODIESEL = sum(exit_PROD.*molar_mass);          % kg/hr
PROD_ANIMFEED = sum(exit_HTLSol.*molar_mass(1:35));   % kg/hr
% ORGANIZE RECYCLE STREAMS
REC_WATER = exit_HTLAq(1)*molar_mass(1);              % kg/hr
REC_NUTRI = REC_WATER*nutri_conc;                     % kg/hr
REC_H2 = sum(exit_H2Rec.*molar_mass);                 % kg/hr
REC_SOLV = C200_Top(35).*molar_mass(35);              % kg/hr
HTL_Oil_Yield = sum(exit_Biocrude.*molar_mass(1:35))/wet_bm(2);



%% PROCESS ENERGY BALANCE MODEL
% [1] ENERGY REQUIREMENTS: Cultivation Process
% (a) Flue Gas Centrifugal Blower Energy (DAYTIME) (in MW) [15]
RATING_FG = 7.5;                          % kW per Blower, 0.12 bar
FG_UNIT_FWRATE = 520;                     % Nm^3/hr capacity
FWRATE_FG = VVM*V_Cult_T*0.06;            % m3/hr Flue Gas Feedrate
NUM_FG = ceil(FWRATE_FG/FG_UNIT_FWRATE);  % Number of Blower Modules
PWR_FG = NUM_FG*RATING_FG;                % kW Power Consumption
ENR_FG = PWR_FG*Daylight*HAR_T*Num_Har*0.0036;  % Daylight FG Bubbling, GJ
% (b) Air Centrifugal Blower Energy (NIGHTTIME) (in MW) [15]
PWR_Air = PWR_FG*(31/73);                 % kW Power Consumption
ENR_Air = PWR_Air*Daylight*HAR_T*Num_Har*0.0036;% Nightime FG Bubbling, GJ
% (c) Stock Solution Pump Energy (in MW) [15]
RATING_Stock = 5.5;                       % kW per Pump, 2.7 bar diffP
FWRATE_Stock = 40;                        % m^3/hr flowrate per pump
NUM_Stock = ceil((MU_Water/(1000*HAR_T))/FWRATE_Stock/STOCK_T);  
PWR_Stock = (NUM_Stock*RATING_Stock)/(NU_Pump*NU_Motor); % kW Power Cons.
ENR_Stock = PWR_Stock*STOCK_T*Num_Har*0.0036; % Stock Transf. Energy, GJ
% (d) Growth Media Pump Energy (in MW) [15]
RATING_GM = 10;                           % kW per Pump, 2.7 bar diffP
FWRATE_GM = 300;                          % m^3/hr flowrate per pump
NUM_GM = ceil((V_Cult_T*Har_Frac/1000)/FWRATE_GM/GM_T);      
PWR_GM = (NUM_GM*RATING_GM)/(NU_Pump*NU_Motor); % kW Power Cons.
ENR_GM = PWR_GM*GM_T*Num_Har*0.0036;      % GM Transf. Energy, GJ 

% [2] ENERGY REQUIREMENTS: Harvest+Dewater Process
% (a) Centrifugal Separator [15]
RATING_SEP = 7.5;                 % kW per Centrifugal Separator
CAPACITY_SEP = 4;                 % Culture processing capacity, m^3/hr
NUM_SEP = ceil((V_Cult_T*Har_Frac/(1000*HAR_T))/CAPACITY_SEP);
PWR_SEP = NUM_SEP*RATING_SEP;             % kW Separation Power Cons.
ENR_SEP = PWR_SEP*SEP_T*Num_Har*0.0036;   % Centrif. separation energy, GJ

% [3] Blowdown Recycle Treatment
% (a) Recycle Growth Media Pump (w/ Filter) [15]
RATING_Rec = 5.5;                 % kW per Recycle Pump
FWRATE_Rec = 40;                  % m^3/hr flowrate processing capacity
NUM_Rec = ceil((recycle_CENTR(1)/1000)/FWRATE_Rec); 
PWR_Rec = (NUM_Rec*RATING_Rec)/(NU_Pump*NU_Motor); % kW Power Consumption
ENR_Rec = PWR_Rec*REC_T*Num_Har*0.0036;   % Broth Recycle pump energy, GJ

% [4] Wet Biomass HTL Treatment
% (a) [ASPEN P-100] Biomass Slurry Compressor Pump
%%%The power consumption for slurry pumps can be calculated by taking the
%%%power consumption for clean water, then multiplying the power by the
%%%Specific Gravity. In addition, the friction losses for clean water must
%%%be multiplied with a correction factor [16] for accurate accounting
FWRATE_Slry = wet_bm(1)/1000;           % m^3/hr slurry flowrate
SG_Slry = SG_MA/(SG_MA - (conc_CENTR*(SG_MA-1)));
PWR_Slry = (FWRATE_Slry*delP_Slry*(1E5/3.6E6)*SG_Slry)/(NU_Pump*NU_Motor);
ENR_Slry = PWR_Slry*HAR_T*Num_Har*0.0036;  % GJ
% (b) [ASPEN H-100/H-200] Slurry PreHeat & HTL Stream Heat Exchanger [20]
REFDUTY_PreHeat = 2800;                 % kW Exchanger Duty w/ Heat Int.
REFRATE_PreHeat = 51360;                % kg/hr Mass Flowrate Balance
PWR_PreHeat = REFDUTY_PreHeat*(sum(comp_biom)/REFRATE_PreHeat);
PWR_PreHeat_Actual = PWR_PreHeat*(1335/2800);
ENR_PreHeat = PWR_PreHeat_Actual*HAR_T*Num_Har*0.0036;  % GJ
% (c) [ASPEN E-100] Hydrothermal Liquefaction Reactor (Black Box)
REFDUTY_HTL = 14.5;                     % kW Reactor Agit.+Pressurization
PWR_HTL = REFDUTY_HTL*(sum(comp_biom)/REFRATE_PreHeat);
ENR_HTL = PWR_HTL*HAR_T*Num_Har*0.0036; % GJ

% [5] Crude Stream Refining Process
% (a) [ASPEN H-300] Crude Oil Phase Pre-Heat for C-100 (Black Box)
REFDUTY_CrudeHeat = 1535.53;            % kW Heater Duty for Crude Pre-Heat
REFRATE_CrudeHeat = 20251.5;            % kg/hr Crude Phase Flowrate
PWR_CrudeHeat = REFDUTY_CrudeHeat*(sum(exit_HTLOil.*molar_mass(1:35))/REFRATE_CrudeHeat);
ENR_CrudeHeat = PWR_CrudeHeat*HAR_T*Num_Har*0.0036;  % GJ
% (b) [ASPEN C-100] Heavy Crude Separation Column (Black Box)
REFDUTY_C100R = 3857.1;                 % kW Reboiler Duty for C-100
REFDUTY_C100C = -3267.3;                % kW Condenser Duty for C-100
REFRATE_C100 = 20251.5;                 % kg/hr Column Flowrate for C-100
PWR_C100R = REFDUTY_C100R*(sum(exit_HTLOil.*molar_mass(1:35))/REFRATE_C100);
PWR_C100C = REFDUTY_C100C*(sum(exit_HTLOil.*molar_mass(1:35))/REFRATE_C100);
ENR_C100R = PWR_C100R*HAR_T*Num_Har*0.0036;  % GJ
ENR_C100C = PWR_C100C*HAR_T*Num_Har*0.0036;  % GJ
% (c) [ASPEN C-200] Light Crude Refining Column (Black Box)
REFDUTY_C200R = 285.7;                  % kW Reboiler Duty for C-200
REFDUTY_C200C = -220.8;                 % kW Condenser Duty for C-200
REFRATE_C200 = 10125.8;                 % kg/hr Column Flowrate for C-200
PWR_C200R = REFDUTY_C200R*(sum(C100_Top.*molar_mass(1:35))/REFRATE_C200);
PWR_C200C = REFDUTY_C200C*(sum(C100_Top.*molar_mass(1:35))/REFRATE_C200);
ENR_C200R = PWR_C200R*HAR_T*Num_Har*0.0036;  % GJ
ENR_C200C = PWR_C200C*HAR_T*Num_Har*0.0036;  % GJ
% (d) [ASPEN C-300] Decane Solvent Recovery Column (Black Box)
REFDUTY_C300R = 898.6;                  % kW Reboiler Duty for C-300
REFDUTY_C300C = -821.9;                 % kW Condenser Duty for C-300
REFRATE_C300 = 9067.6;                  % kg/hr Column Flowrate for C-300
PWR_C300R = REFDUTY_C300R*(sum(C200_Bot.*molar_mass(1:35))/REFRATE_C300);
PWR_C300C = REFDUTY_C300C*(sum(C200_Bot.*molar_mass(1:35))/REFRATE_C300);
ENR_C300R = PWR_C300R*HAR_T*Num_Har*0.0036;  % GJ
ENR_C300C = PWR_C300C*HAR_T*Num_Har*0.0036;  % GJ

% [6] Oil Stream Hydrotreating Process
% (a) [ASPEN P-300/H-400] Biocrude Compressor Pump
REFDUTY_Upg = 1995.5;                   % kW Heat + Compressor Duty 
REFRATE_Upg = 14160.4;                  % kg/hr Biocrude Stream
PWR_Upg = REFDUTY_Upg*(sum(exit_Biocrude.*molar_mass(1:35))/REFRATE_Upg);
ENR_Upg = PWR_Upg*HAR_T*Num_Har*0.0036;
% (b) [ASPEN P-400] Hydrogen Compressor Pump
FWRATE_H2 = UPG_In(14)*24;              % m^3/hr Hydrogen feedrate
delP_H2 = 7.9;                          % Pressure change (bar) for Upgrading
PWR_H2 = (FWRATE_H2*delP_H2*(1E5/3.6E6))/(NU_Pump*NU_Motor);
ENR_H2 = PWR_H2*HAR_T*Num_Har*0.0036;   % GJ
% (c) [ASPEN F-100] Phase Separating Decanter
REFDUTY_Fract = -5.53;                  % kW Decanter Duty for Temp Control
REFRATE_Fract = 14260.4;                % kg/hr Upgraded Biocrude Flowrate
PWR_Fract = REFDUTY_Fract*(sum(UPG_Out.*molar_mass)/REFRATE_Fract);
ENR_Fract = PWR_Fract*HAR_T*Num_Har*0.0036; % GJ

% ORGANIZE ENERGY DUTIES FOR UTILITY CALCULATIONS
DUTY_Elec = ENR_FG+ENR_Air+ENR_Stock+ENR_GM+ENR_SEP+ENR_Rec+ENR_Slry+ENR_HTL+ENR_Upg+ENR_H2;
DUTY_HPSteam = ENR_PreHeat;
DUTY_MPSteam = ENR_CrudeHeat+ENR_C100R+ENR_C200R+ENR_C300R;
DUTY_CoolW = ENR_C100C+ENR_C200C+ENR_C300C+ENR_Fract;



%% TECHNOECONOMIC EVALUATION [TEA] MODEL
%============================ CAPITAL EXPENSES ===========================%
% [1] Cultivation Process
% (a) Vertical Airlift Plastic Bag Cultivators [21]
%%%One unit consists of 24 panels each with volume 48m*0.7m*0.045m at fill
%%%Based on Algenol's low cost hanging bag airlift PBRS reported by [21]
%%%Each cultivation unit of 24 panels occupies 1250m^2 of land area
Num_Cult_Unit = ceil((V_Cult_T/1000)/V_Cult_Unit);    % Number of PBR Units
T_Land = Land_Cult_Unit*Num_Cult_Unit;                % Total Land Occupied
EQ_PBR_Cult = C_PBR*(T_Land/4046.86)*(CEPCI/CEPCI_PBR);
% (b) Flue Gas Centrifugal Blowers
EQ_FG_Blower = C_FG_Blower*NUM_FG*(CEPCI/CEPCI_PumpsBlowers);    
% (c) Stock Solution Centrifugal Pumps
EQ_Stock_Pump = C_Stock_Pump*NUM_Stock*(CEPCI/CEPCI_PumpsBlowers);
% (d) Growth Media Centrifugal Pumps
EQ_GM_Pump = C_GM_Pump*NUM_GM*(CEPCI/CEPCI_PumpsBlowers);       

% [2] Harvest+Dewater Process
% (a) Centrifugal Heavy Duty Separator 
EQ_SEP = C_SEP*NUM_SEP*(CEPCI/CEPCI_Separator)/SEP_T;

% [3] Blowdown + Media Recycle
% (a) Recycle Growth Media Centrifugal Pump
EQ_Rec_Pump = C_Rec_Pump*NUM_Rec*(CEPCI/CEPCI_PumpsBlowers);
% (b) Stock Solution Holding Tanks (Small Field Erected)
NUM_SS_Tank = ceil((MU_Water/1000)/CAP_SS_T);
EQ_SS_Tank = REFCOST_SS_T*NUM_SS_Tank*(CEPCI/CEPCI_Upstr_Tanks);
% (c) Growth Media Holding Tanks (Large Field Erected)
NUM_GM_Tank = ceil((V_Cult_T*Har_Frac/1000)/CAP_GM_T);
EQ_GM_Tank = REFCOST_GM_T*NUM_GM_Tank*(CEPCI/CEPCI_Upstr_Tanks);

% [4] Wet Biomass HTL Treatment
% (a) [ASPEN P-100] Biomass Slurry Compressor Pump (Nth Plant)
EQ_P100 = REFCOST_P100*(CEPCI/CEPCI_DS)*((wet_bm(1)/1000)/116.75)^NTH_PMP;
% (b) [ASPEN H-100/H-200] Slurry PreHeat & HTL Stream Heat Exchanger [20]
EQ_H100 = REFCOST_H100*(CEPCI/CEPCI_DS)*(sum(comp_biom)/73155.4)^NTH_HEX;
% (c) [ASPEN E-100] Hydrothermal Liquefaction Reactor (Nth Plant)
EQ_E100 = REFCOST_E100*(CEPCI/CEPCI_DS)*(wet_bm(2)/25333.3)^NTH_CSTR;

% [5] Crude Stream Refining Process
% (a) [ASPEN H-300] Crude Oil Phase Pre-Heat for C-100 (Nth Plant)
EQ_H300 = REFCOST_H300*(CEPCI/CEPCI_DS)*(sum(exit_HTLOil.*molar_mass(1:35))/REFRATE_CrudeHeat)^NTH_DIST;
% (b) [ASPEN C-100] Heavy Crude Separation Column (Nth Plant)
EQ_C100 = REFCOST_C100*(CEPCI/CEPCI_DS)*(sum(exit_HTLOil.*molar_mass(1:35))/REFRATE_C100)^NTH_DIST;
% (c) [ASPEN C-200] Light Crude Refining Column (Nth Plant)
EQ_C200 = REFCOST_C200*(CEPCI/CEPCI_DS)*(sum(C100_Top.*molar_mass(1:35))/REFRATE_C200)^NTH_DIST;
% (d) [ASPEN C-300] Decane Solvent Recovery Column (Nth Plant)
EQ_C300 = REFCOST_C300*(CEPCI/CEPCI_DS)*(sum(C200_Bot.*molar_mass(1:35))/REFRATE_C300)^NTH_DIST;

% [6] Oil Stream Hydrotreating Process
% (a) [ASPEN P-300/H-400] Biocrude Compressor Pump
EQ_P300 = REFCOST_P300*(CEPCI/CEPCI_DS)*(sum(exit_Biocrude.*molar_mass(1:35))/REFRATE_Upg)^NTH_PMP;
% (b) [ASPEN P-400] Hydrogen Compressor Pump
EQ_P400 = REFCOST_P400*(CEPCI/CEPCI_DS)*(UPG_In(14)*24/1190.55)^NTH_COMP;
% (c) [ASPEN F-100] Phase Separating Decanter
EQ_F100 = REFCOST_F100*(CEPCI/CEPCI_DS)*(sum(UPG_Out.*molar_mass)/REFRATE_Fract)^NTH_DIST;

% CAPITAL COST SUMMARY
ISBL_Cult = EQ_PBR_Cult+EQ_FG_Blower+EQ_Stock_Pump+EQ_GM_Pump;
ISBL_HarDew = EQ_SEP+EQ_Rec_Pump+EQ_SS_Tank+EQ_GM_Tank;
ISBL_HTL = EQ_P100+EQ_H100+EQ_E100+EQ_H300+EQ_C100+EQ_C200+EQ_C300;
ISBL_UPG = EQ_P300 + EQ_P400 + EQ_F100;
ISBL_Total = ISBL_Cult + ISBL_HarDew + ISBL_HTL + ISBL_UPG;
FIXED_CAPEX = ISBL_Total*(1+OSBL_OS)*(1+OSBL_DE+OSBL_CN);
CRF = (DISCO*(DISCO+1)^LIFET)/(((DISCO+1)^LIFET)-1);
ANNUALIZED_CAPEX = CRF*FIXED_CAPEX;

%=========================== OPERATING EXPENSES ==========================%
% [1] Plant Expenses
% (a) Land, Clean-In-Place & PBR Plastic Replacement Costs
COST_PE_LCP = (OP_PE_LAND+OP_PE_CIP+OP_PE_PBR)*(T_Land/4046.86);
% (b) Innoculum Costs
COST_PE_INOC = sum([OP_RM_NSOURCE, OP_RM_PSOURCE, OP_RM_SSOURCE].*(nutri_conc*V_Cult_T))/LIFET;

% [2] Raw Material Costs
% (a) Flue Gas Utilization Costs
MOLES_FG = (prod_bm*cult_stoic(2)*Num_Har*(1/44.01))/P_CO2;
COST_RM_FG = OP_RM_FG*MOLES_FG*sum(FG_MM.*FG_Moles);
% (b) Nutrient Stock Make-Up Costs
MU_Nutri = (MU_Nutri-(REC_NUTRI*HAR_T))*Num_Har;  % kg/yr Makeup Nutrients
COST_RM_NU = sum([OP_RM_NSOURCE, OP_RM_PSOURCE, OP_RM_SSOURCE].*MU_Nutri);
% (c) Water Make-up Costs 
MU_Water = (MU_Water-(REC_WATER*HAR_T))*Num_Har;  % kg/yr Makeup Water
COST_RM_H2O = OP_RM_H2O*MU_Water;             
% (d) Oil-Phase Solvent Make-Up Costs
MU_Solv = (ExSolv*molar_mass(35)-REC_SOLV)*HAR_T*Num_Har;
COST_RM_SOLV = OP_RM_SOLV*(MU_Solv + (ExSolv*molar_mass(35))/LIFET);
% (e) Hydrogen Make-Up Costs
MU_H2 = (UPG_In(14)*molar_mass(14)-REC_H2)*HAR_T*Num_Har;
COST_RM_H2 = OP_RM_H2*(MU_H2 + (UPG_In(14)*molar_mass(14))/LIFET);

% [3] Utility Costs
COST_UT_ELEC = UT_ELEC*DUTY_Elec;          % Yearly Electrical Costs
COST_UT_HPS = UT_HP_STEAM*DUTY_HPSteam;    % Yearly HP Steam Costs
COST_UT_MPS = UT_MP_STEAM*DUTY_MPSteam;    % Yearly MP Steam Costs
COST_UT_CW = UT_CW*DUTY_CoolW;             % Yearly Chilled Water Costs

% OPERATING COST SUMMARY
TOTAL_OPEX_PE = COST_PE_LCP + COST_PE_INOC;
TOTAL_OPEX_RM = COST_RM_FG + COST_RM_NU + COST_RM_H2O + COST_RM_SOLV + COST_RM_H2;
TOTAL_OPEX_UT = COST_UT_ELEC + COST_UT_HPS + COST_UT_MPS + COST_UT_CW;

%============================== COST SUMMARY =============================%
ANNUAL_PROD_COST = ANNUALIZED_CAPEX + TOTAL_OPEX_PE + TOTAL_OPEX_RM + TOTAL_OPEX_UT;
BYPROD_ANIMFEED = 1.85;                    % Avg. $/kg sell price from [26]
SPECIFIC_PROD_COST = ANNUAL_PROD_COST/(PROD_BIODIESEL*HAR_T*Num_Har);% $/kg
SPECIFIC_CREDITS = PROD_ANIMFEED*BYPROD_ANIMFEED/PROD_BIODIESEL;
NET_SPECIFIC_PROD_COST = SPECIFIC_PROD_COST - SPECIFIC_CREDITS;



%% CO2 LIFE CYCLE ASSESSMENT MODEL
% [1] Direct Plant Emissions of GHG
% (a) Cultivation Off Gases **NOT COUNTED >> ASSUME GO TO CAPTURE PLANT**
CULT_CO2_In_Moles = (VVM*V_Cult_T*0.06/RTP_FG)*(P_CO2/1000)*HAR_T*Num_Har;
CULT_CO2_In = CULT_CO2_In_Moles*FG_MM(1);      % kg CO2 in Total
CULT_CO2_Consumed = cons_cult(2)*Num_Har;
% (b) Sequestration as TIC in Media
MOLES_CO2_SEQ = cult_out(NUM_INT,2)*V_Cult_T;  % CO2 sequestered as TIC
MASS_CO2_SEQ = MOLES_CO2_SEQ*(FG_MM(1)/1000)*Num_Har; % Total Sequestered
CULT_CO2_Out = CULT_CO2_In - CULT_CO2_Consumed - MASS_CO2_SEQ;
% (c) HTL Reactor Gas Phase Emissions
HTL_CO2_Out = exit_HTLGas(13)*molar_mass(13)*HAR_T*Num_Har;
% Sum all Direct Emissions
%TOT_EM_DIR = CULT_CO2_Out + HTL_CO2_Out;          % MASS BALANCE APPROACH
TOT_EM_DIR = -CULT_CO2_Consumed + HTL_CO2_Out;     % MITIGATION APPROACH

% [2] Indirect Plant Emissions of GHG
% (a) Indirect Emissions from Energy Consumption
EM_ENR_ELEC = GWI_ELEC_GRID*DUTY_Elec;
EM_ENR_STEAM_HP = GWI_STEAM_HP*DUTY_HPSteam;
EM_ENR_STEAM_MP = GWI_STEAM_MP*DUTY_MPSteam;
EM_ENR_COOLW = GWI_COOLW*DUTY_CoolW;
TOT_EM_ENR = EM_ENR_ELEC + EM_ENR_STEAM_HP + EM_ENR_STEAM_MP + EM_ENR_COOLW;
% (b) Indirect Emissions from Raw Material Consumption
EM_RM_NSOURCE = (((nutri_conc(1)*V_Cult_T)/LIFET)+MU_Nutri(1))*GWI_NSOURCE;
EM_RM_PSOURCE = (((nutri_conc(2)*V_Cult_T)/LIFET)+MU_Nutri(2))*GWI_PSOURCE;
EM_RM_SSOURCE = (((nutri_conc(3)*V_Cult_T)/LIFET)+MU_Nutri(3))*GWI_SSOURCE;
EM_RM_WATER = MU_Water*GWI_WATER;
EM_RM_SOLV = (MU_Solv + (ExSolv*molar_mass(35))/LIFET)*GWI_SOLV;
EM_RM_H2 = (MU_H2 + (UPG_In(14)*molar_mass(14))/LIFET)*GWI_H2;
TOT_EM_RM = EM_RM_NSOURCE + EM_RM_PSOURCE + EM_RM_SSOURCE + EM_RM_WATER + EM_RM_SOLV + EM_RM_H2;
% (c) Indirect Emissions from Plant Construction and Salvage
TOT_EM_CaS = 0;             % Not Studied
% (d) Indirect Emissions from Product Consumption
TOT_EM_PC = 0;              % Cradle to Gate (Product Production) only

% CO2 GWI SUMMARY (in kg CO2-eq)
TOTAL_GWI = B_DIR*TOT_EM_DIR + B_ENR*TOT_EM_ENR + B_MAT*TOT_EM_RM + B_CaS*TOT_EM_CaS + B_PC*TOT_EM_PC;
SPECIFIC_GWI = TOTAL_GWI/(PROD_BIODIESEL*HAR_T*Num_Har);



%% MASS BALANCE INFO FOR EXPORT %%
FG_Moles_Total = (VVM*V_Cult_T*0.06/RTP_FG)*(P_CO2/1000)/P_CO2;
FG_Moles_Flowrate = FG_Moles_Total.*FG_Moles;
Stream1and2_Flowrate = FG_Moles_Flowrate.*FG_MM;   
Stream3_Flowrate = MU_Water/(HAR_T*Num_Har);
Stream4_Flowrate = MU_Nutri/(HAR_T*Num_Har);
Stream10_Flowrate = blowdown_loss;   % Water, Biomass, N, P, S
Stream11_Flowrate = cult_recycle;
Stream12_Flowrate = wet_bm;
Stream7_Flowrate = [Stream1and2_Flowrate(1) - (MASS_CO2_SEQ/(HAR_T*Num_Har)), Stream1and2_Flowrate(2:5), Stream1and2_Flowrate(6) + (prod_bm*2.05/HAR_T), Stream1and2_Flowrate(7)];
CO2_Fixation_Rate = CULT_CO2_Consumed/(Num_Har*HAR_T);
Stream8_Flowrate = exit_cult;



%% EXPORT EVALUATION METRICS
EvalMetrics = [NET_SPECIFIC_PROD_COST, SPECIFIC_GWI];      
end