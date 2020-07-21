%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% PARAMETER RANKING via MC-NAIVE BAYES %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Jeehwan Lee, KAIST, stevelee@kaist.ac.kr %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clc;
tic;
format longG;
progress = waitbar(0, 'Running...', 'Name', 'Running PR-MC-Sobol...');
total_steps = 1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [1] Initialize Parameters, Spaces, and Bayesian Targets
% (a) Define # of Model Evaluation (# of MC sampling of Parameter Sets)
n_sim = 2500;
export = 1;
if n_sim > 50000
    simwarn = questdlg('WARNING: For n_sim >50,000 there could be issues with exporting results via xlswrite. Disable xlswrite?',...
        'WARNING',...
        'No','Yes','Yes');
    switch simwarn
        case 'Yes'
            export = 0;
        case 'No'
            export = 1;
    end
end

% (b) Automatically generate Bayesian target values, or manually enter them
bayes_auto = 1;
bayes_target = questdlg('Automatically generate Target CCU Output Metric Values for Bayesian Kernel Distribution?',...
        'USER INPUT',...
        'No','Yes','Yes');
switch bayes_target
    case 'Yes'
        bayes_auto = 1;
    case 'No'
        bayes_auto = 0;
end


% (c) Initialize Parameter Values for Key Parameters with kp
%     The number of elements in kp should correspond to # of key parameters 
%     that are considered in the model. If the # of key params of the model 
%     are changed, the length of kp should be adjusted accordingly.
%     GUI Support will be added in the future
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kp = [1 1 1 1 1 1 1 1 1 1];                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           

% (d) Scalar # of Key Parameters. Count columns of kp array.
n_kp = size(kp, 2);

% (e) Define two identical parameter spaces for Sobol Analysis
%     If n_kp params are simulated n_sim times, space is n_sim x n_kp 
ParSpace = [];                 % Space for the Key Parameter of Interest

% (f) Define Resolution of Kernel Distributions
np_1 = 10000;               % Resolution for Parameter KDs for MC sampling
np_2 = 50000;               % Resolution for Bayesian KDs for computing area difference

%%Progress Bar%%
waitbar(20/total_steps, progress, 'Generating Parameter Spaces...');



%% [2] Populate Parameter Space via MC Simulation
% Distributions can either be defined parametrically or non-parametrically. 
% If  repeated experiment sampling under identical conditions is possible, 
% then a parametric sample distribution can be defined with mean and stdev.
% Otherwise, a PDF can be generated non-parametrically via kernel density
% [EX] Parametric Distribution:
% ParSpace(:,1) = 3.4.*randn(n_sim, 1) + 0.9                 
% ParSpace(:,3) = 10.*randn(n_sim, 1) + 24.1                
%
% [EX] Kernel Distribution using fitdist
% data_x1 = [13.5, 16.4, 14.4, 19.6, 15.0, 15.9, 16.2]  % Sample data
%
% [EX] Kernel Distribution using ksdensity (slower than fitdist!)
% ParSpace(:,1) = 35.*randn(n_sim, 1) + 325;     % Parametric PDF for x1
% data_x2 = [13.5, 16.4, 14.4, 19.6, 15.0, 15.9, 16.2];
% data_x3 = [82.4, 77.6, 83.3, 80.1];
% [f_x2, xi_x2] = ksdensity(data_x2, 'npoints', np_1);
% [f_x3, xi_x3] = ksdensity(data_x3, 'npoints', np_1);
% for i = 1:n_sim                                  % Kernel PDFs for x2,x3
%     ParSpace(i,2) = randarb(xi_x2, f_x2);
%     ParSpace(i,3) = randarb(xi_x3, f_x3);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Define Parametric Distributions for Key Params]
ParSpace(:,1) = 0.0016.*randn(n_sim, 1) + 0.0224;               % avg_GPC
ParSpace(:,8) = 13030.6.*randn(n_sim, 1) + 65153;               % C_PBR

% [List Kernel Data Points]
data_x2 = [0.89, 1.07, 1.24, 1.00, 0.4, 0.8, 0.84,...
           0.81, 1.11, 0.82, 1];           % Kernel Data for Mu_max
data_x3 = [3.6, 2.042, 1.4];               % Kernel Data for kLa
data_x4 = [0.3, 0.14, 0.18, 0.22, 0.1218,...
           0.19, 0.201, 0.225];            % Kernel Data for y_Lipid
data_x5 = [0.1358, 0.1184, 0.1283, 0.116]; % Kernel Data for P_CO2
data_x6 = [0.402, 0.424, 0.465, 0.515, 0.562, 0.592, 0.599, 0.581,...
           0.542, 0.493, 0.445, 0.411];    % Kernel Data for Daylight %
data_x7 = [0.95, 0.93, 0.99, 0.97];        % Kernel Data for Recovery %
data_x9 = [187.5, 247.2, 278.06, 313.88, 469.17, 313.06, 288.88, 250,...
           266.11, 241.66];                % Kernel Data for Elec. GWI
data_x10 = [3.234, 4.62, 6.006, 3.469, 3.528, 3.527, 6.2, 2.74, 1.59,...
            1.13];                         % Kernel Data for N-Fert. GWI

% [Define Kernel Distributions]
f_x2 = fitdist(data_x2', 'Kernel', 'Support', [0.01, 10]);
f_x3 = fitdist(data_x3', 'Kernel', 'Support', [0.01, 100]);
f_x4 = fitdist(data_x4', 'Kernel', 'Support', [0.01, 0.9]);
f_x5 = fitdist(data_x5', 'Kernel', 'Support', [0.01, 0.9]);
f_x6 = fitdist(data_x6', 'Kernel', 'Support', [0.25, 0.75]);
f_x7 = fitdist(data_x7', 'Kernel', 'Support', [0.7, 1]);
f_x9 = fitdist(data_x9', 'Kernel', 'Support', [1, 999]);
f_x10 = fitdist(data_x10', 'Kernel', 'Support', [0.01, 100]);

% [Sample from Kernel Distributions]
% Random Sample to Populate Parameter Space
ParSpace(:,2) = random(f_x2, n_sim, 1);
ParSpace(:,3) = random(f_x3, n_sim, 1);
ParSpace(:,4) = random(f_x4, n_sim, 1);
ParSpace(:,5) = random(f_x5, n_sim, 1);
ParSpace(:,6) = random(f_x6, n_sim, 1);
ParSpace(:,7) = random(f_x7, n_sim, 1);
ParSpace(:,9) = random(f_x9, n_sim, 1);
ParSpace(:,10) = random(f_x10, n_sim, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Progress Bar%%
waitbar(100/total_steps, progress, 'Calculating Monte Carlo Outputs');



%% [3] Process Parameter Space for CCU Evaluation
% (a) Determine the # of CCU Model Outputs (# of Evaluation Metrics)
n_out = size(CCUS_Biocrude(kp),2);

% (b) Evaluate CCU Model from MC Parameter Space
%     For n_out outputs, we need n_sim*n_out evaluations of CCU model
outputs = zeros(n_sim, n_out);
parfor i = 1:n_sim
    % Each Parameter Set is a Row in the ParamSpace matrix
    Parameter_Set = ParSpace(i,:);   
    for j = 1:n_out
        output_temp = CCUS_Biocrude(Parameter_Set);
        outputs(i,j) = output_temp(j);
    end
end
%%Progress Bar%%
waitbar(250/total_steps, progress, 'Calculating Monte Carlo Outputs');



%% [4] Bayesian Classification of CCU Model Outputs
% (a) Define Target Values for each CCU Evaluation Metric
% [EX] If the CCU model has 3 output metrics, i.e., n_out=3 then "target" 
% must be an array of length 3 with the target values for those metrics. 
% NOTE: If the target values are too low/high compared to that of MC 
% simulated outputs the Bayesian difference for "Success" and "Failure" 
% becomes large, which might make the Bayesian classification process 
% unreliable. It is best to choose realistic target values so that there is
% some overlap between the kernel distributions
% NOTE: Must be an array of n_out length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
targets = zeros(1, n_out);  % Array of Target Bayesian Output Metric Values
if bayes_auto == 1
    for i = 1:n_out
        targets(i) = mean(outputs(:,i));
    end
else
    for i = 1:n_out
        prompt = sprintf('Enter Target Value for Model Output Metric %.0f \n', i);
        prompt_title = 'Enter Target Values';
        dims = [1 35];
        answer = inputdlg(prompt, prompt_title, dims);
        targets(i) = str2double(answer{1});
    end
end
%targets = [4.92, 3.33];  % Unit Prod Cost($/kg), Specific GWI(kg CO2eq/kg)              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (b) Generate Bayesian Classification Matrix
% A MC-sampled output is classified as a "success" if the value of the 
% output metric is less than or equal to the target metric value above. An 
% output is classified as a failure if the evaluation metric value has a 
% higher value than the target. 
% NOTE: The signs must be consistent across all metrics. We propose using 
% cost metrics such as CAPEX, Unit Produciton Cost, and Specific GWI which 
% is better the smaller the value. If the CCU model includes profitability 
% metrics (ex: Revenue, Net Profit), the code below must be adjusted.
classmat = zeros(n_sim, n_out);
parfor i = 1:n_sim
    for j = 1:n_out
        if outputs(i,j) <= targets(j)
            classmat(i,j) = 1;
        else
            classmat(i,j) = 0;
        end
    end
end
%%Progress Bar%%
waitbar(500/total_steps, progress, 'Generating Bayesian Classification Matrix');



%% [5] Factorize Classification Matrix and Append Parameter Values
% (a) Generate Empty Factorized Containers. The maximum dimensions of each
%     Label matrix is n_sim x n_kp x n_out
Success_Mat = zeros(n_sim, n_kp, n_out);
Failure_Mat = zeros(n_sim, n_kp, n_out);

% (b) Loop through "classmat" then append values from ParamSpace to 
%     Success_Mat or Failure_Mat accordingly
parfor i = 1:n_sim
    for j = 1:n_out
        if classmat(i,j) == 1
            Success_Mat(i,:,j) = ParSpace(i,:);
            Failure_Mat(i,:,j) = zeros(1,n_kp,'uint32');
        else
            Success_Mat(i,:,j) = zeros(1,n_kp,'uint32');
            Failure_Mat(i,:,j) = ParSpace(i,:);
        end
    end
end
%%Progress Bar%%
waitbar(700/total_steps, progress, 'Generating Bayesian Kernel Distributions');



%% [6] For each KP in each CCU Eval Metric, generate Kernel Distributions
% (a) Initialize Containers
S_KernelData = zeros(n_sim, n_kp, n_out); % Sampled matrix for "success" 
F_KernelData = zeros(n_sim, n_kp, n_out); % Sampled matrix for "failure" 
NBC_KDE_Rank = zeros(n_kp, n_out);        % A n_kp x n_out matrix of ranks


% (b) Populate NBC_KDE_Rank matrix by evaluating the differences
parfor i = 1:n_kp
    for j = 1:n_out
        % Generate array of sample data for Bayesian Kernel Distributions
        S_KernelData = nonzeros(Success_Mat(:,i,j))';    
        F_KernelData = nonzeros(Failure_Mat(:,i,j))';
        % Generate the Kernel Distribution for each Bayesian Outcome
        [F_s, x_s] = ksdensity(S_KernelData, 'npoints', np_2);
        [F_f, x_f] = ksdensity(F_KernelData, 'npoints', np_2);
        % Calculate PDF difference between "Success" & "Failure" Kernels
        NBC_KDE_Rank(i,j)  = trapz(x_s, abs(F_s-F_f));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Plot the Kernel Distributions (Less overlap = Smaller Score)
        %figure(i)
        %plot(x_s,F_s,'LineWidth',2);
        %hold on
        %plot(x_f,F_f,'r--','LineWidth',2);
        %title(['For Parameter ', num2str(i), ' in output ', num2str(j)])
        %hold off
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%%Progress Bar%%
waitbar(950/total_steps, progress, 'Sorting Parameters by NBC-KDE Rank');



%% [7] Rank Parameters and Export Data
% (a) Rank Parameters by sorting the NBC-KDE matrix
parfor i = 1:n_out
    [score, rank] = sort(NBC_KDE_Rank(:,i), 'descend');
    fprintf('Parameters Ordered by NBC Rank for CCU Eval Output %.0f \n',i)
    rank
    sprintf('The Rank Scores for the Above Order are: ')
    NBC_KDE_Rank(:,i)
end

% (b) Generate plots (Example for via For loop below)
nbins = 30;             % Define resolution of output metric's historgram
%for i = 1:n_out
    %figure(i)
    %histogram(outputs(:,i), nbins, 'facecolor', [1/i, 1/(4*i), 1/(16*i)])
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Output #1
    figure('Name', 'Variance of Unit Prod Cost')
    histogram(outputs(:,1), nbins, 'facecolor', [0, 0, 0])
    xlabel('Biocrude Unit Prod. Cost, $/kg')
    ylabel('Frequency')
    % Output #2
    figure('Name', 'Variance of Specific GWI')
    histogram(outputs(:,2), nbins, 'facecolor', [0, 0.5, 0])
    xlabel('Specific GWI, kg-CO2-eq/kg')
    ylabel('Frequency')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
targets

% (c) Export Results
%%Progress Bar%%
waitbar(1000/total_steps, progress, 'Complete');
delete(progress)
toc

if export == 1
    xlswrite('MC-NBC_ParameterSpace.xls', ParSpace)
    xlswrite('MC-NBC_EvaluatedOutputs.xls', outputs)
    xlswrite('MC-NBC_EvalMetricTargets.xls', targets)
    xlswrite('MC-NBC_ClassifiedOutputs.xls', classmat)
    xlswrite('MC-NBC_SuccessKernelData.xls', S_KernelData)
    xlswrite('MC-NBC_FailureKernelData.xls', F_KernelData)
    xlswrite('MC-NBC_Parameters_Ranked.xls', NBC_KDE_Rank)
end