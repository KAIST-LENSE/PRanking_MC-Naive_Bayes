% This is just a trial function for understanding Parameter Ranking Methods
function EvalMetrics = SampleModel(kp)             %kp = key parameters
%% Key Process Parameters
% # of Params should equal CCUS model DOF
% Key Params are the independent variables of the process model, such as Reflux Ratio,
% bottoms rate, feed flowrate, CAPEX factors, etc
x1 = kp(1);                 % 1st key parameter
x2 = kp(2);                 % 2nd key parameter
x3 = kp(3);                 % 3rd key parameter
%...etc

%% Process Model
Yield = 2*x1 + 2*x2;
Energy_Elec = 3*x1;
Energy_Steam = x2 + x3;
Size_Reactor = x1-x3;
Size_Column = 4*(x1/x2);
%...etc

%% Model Outputs
CAPEX = 2*Size_Reactor + Size_Column;
OPEX = Energy_Elec + Energy_Steam;
Unit_Prod_Cost = CAPEX + OPEX;
Specific_GWI = 3*Energy_Elec + Energy_Steam;

%% Output Array
EvalMetrics = [Unit_Prod_Cost, Specific_GWI];      
end