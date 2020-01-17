%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Parameter Ranking via MC-Naive Bayes %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ***REQUIRES MATLAB STATISTICS TOOLBOX!!*** %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Jeehwan Lee, KAIST, stevelee@kaist.ac.kr %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;
tic;
%% [1] Initialize Parameters and Spaces
% (a) Define # of Model Evaluation (# of MC sampling of Parameter Sets)
n_sim = 5000;

% (b) Initialize Parameter Values for Key Parameters with kp
%     The number of elements in kp should correspond to # of key parameters that are
%     considered in the CCU model. If the # of key params of the CCU model are 
%     changed, the length of kp should be adjusted accordingly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kp = [1 1 1];                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (c) Scalar # of Key Parameters. Count columns of kp array.
n_kp = size(kp, 2);

% (d) Define the Parameter Space
%     If n_kp params are simulated n_sim times, a space of n_sim x n_kp is needed
ParamSpace = [];       

% (e) Define Resolution of Kernel Distributions
np_1 = 10000;               % Resolution for Parameter KDs for MC sampling
np_2 = 50000;               % Resolution for Bayesian KDs for computing area difference


%% [2] Populate Parameter Space via MC Simulation
% Distributions can either be defined parametrically or non-parametrically. If 
% repeated experiment sampling under identical conditions is possible, then a  
% parametric sample distribution can be defined with mean and stdev. Otherwise, a PDF
% can be generated non-parametrically using sample data + Gaussian kernels. See below
% [EX] Parametric Distribution:
% ParamSpace(:,1) = 3.4.*randn(n_sim, 1) + 0.9                 >>> Mean=0.9, Stdev=3.4
% ParamSpace(:,3) = 10.*randn(n_sim, 1) + 24.1                 >>> Mean=24.1, Stdev=10
%
% [EX] Non-Parametric (Kernel) Distribution:
% data_var3 = [13.5, 16.4, 14.4, 19.6, 15.0, 15.9, 16.2]..etc  >>> Sample data points
% [f_var3, xi_var3] = ksdensity(data_var3, 'npoints', np_1)    >>> Generate KDE PDF
% for i = 1:n_sim                                              >>> Sample from KDE PDF
%     ParamSpace(i,3) = randarb(xi_var3, f_var3);                  using "randarb"
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ParamSpace(:,1) = 35.*randn(n_sim, 1) + 325;                  % Parametric PDF for x1
data_x2 = [13.5, 16.4, 14.4, 19.6, 15.0, 15.9, 16.2];
data_x3 = [82.4, 77.6, 83.3, 80.1];
[f_x2, xi_x2] = ksdensity(data_x2, 'npoints', np_1);
[f_x3, xi_x3] = ksdensity(data_x3, 'npoints', np_1);
for i = 1:n_sim                                               % Kernel PDFs for x2, x3
    ParamSpace(i,2) = randarb(xi_x2, f_x2);
    ParamSpace(i,3) = randarb(xi_x3, f_x3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [3] Evaluate CCU Model using MC sampled Parameter Space
% (a) Determine the # of CCU Model Outputs (# of Evaluation Metrics)
n_out = size(SampleModel(kp),2);

% (b) Evaluate CCU Model from MC Parameter Space
%     For n_out outputs, we need n_sim*n_out evaluations of CCU model
outputs = [];
for i = 1:n_sim
    % Each Parameter Set is a Row in the ParamSpace matrix
    Parameter_Set = ParamSpace(i,:);   
    for j = 1:n_out
        output_temp = SampleModel(Parameter_Set);
        outputs(i,j) = output_temp(j);
    end
end


%% [4] Bayesian Classification of CCU Model Outputs
% (a) Define Target Values for each CCU Evaluation Metric
% [EX] If the CCU model has 3 output metrics, i.e., n_out=3 (Unit Production Cost,
%      CAPEX, Specific GWI) then "target" must be an array of length 3 with the target
%      values for those metrics. 
% NOTE: If the target values are too low/high compared to that of MC simulated outputs
% the Bayesian difference for "Success" and "Failure" becomes large, which might make
% the Bayesian classification process unreliable. It is best to choose realistic
% target values so that there is some overlap between the kernel distributions
% NOTE: Must be an array of n_out length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
targets = [360, 32];                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (b) Generate Bayesian Classification Matrix
% A MC-sampled output is classified as a "success" if the value of the output metric
% (ex: CAPEX) is less than or equal to the target metric value above. An output is
% classified as a failure if the evaluation metric value has a higher value than the
% target. 
% NOTE: The signs must be consistent across all metrics. We propose using cost metrics
% such as CAPEX, Unit Produciton Cost, and Specific GWI which is better the smaller
% the value. If the CCU model includes profitability metrics (ex: Revenue, Net
% Profit), the code below must be adjusted accordingly
classmat = [];
for i = 1:n_sim
    for j = 1:n_out
        if outputs(i,j) <= targets(j);
            classmat(i,j) = 1;
        else
            classmat(i,j) = 0;
        end
    end
end


%% [5] Factorize Classification Matrix and Append Parameter Values
% (a) Generate Empty Factorized Containers. The maximum dimensions of each Label
%     matrix is n_sim x n_kp x n_out
Success_Mat = [];
Failure_Mat = [];

% (b) Loop through "classmat" then append values from ParamSpace to Success_Mat or
%     Failure_Mat accordingly
for i = 1:n_sim
    for j = 1:n_out
        if classmat(i,j) == 1
            Success_Mat(i,:,j) = ParamSpace(i,:);
            Failure_Mat(i,:,j) = zeros(1,n_kp,'uint32');
        else
            Success_Mat(i,:,j) = zeros(1,n_kp,'uint32');
            Failure_Mat(i,:,j) = ParamSpace(i,:);
        end
    end
end


%% [6] For each KP in each CCU Eval Metric, generate Kernel Distributions
% (a) Initialize Containers
S_KernelData = [];      % Temporary matrix with sample data for "success" kernel distr
F_KernelData = [];      % Temporary matrix with sample data for "failure" kernel distr
NBC_KDE_Rank = [];      % A n_kp x n_out matrix of column rank data

% (b) Populate NBC_KDE_Rank matrix by evaluating the differences in the Bayesian PDFs
for i = 1:n_kp
    for j = 1:n_out
        % Generate array of sample data for Bayesian Kernel Distributions
        S_KernelData = nonzeros(Success_Mat(:,i,j))';    
        F_KernelData = nonzeros(Failure_Mat(:,i,j))';
        % Generate the Kernel Distribution for each Bayesian Outcome
        [F_s, x_s] = ksdensity(S_KernelData, 'npoints', np_2);
        [F_f, x_f] = ksdensity(F_KernelData, 'npoints', np_2);
        % Calculate the PDF difference between "Success" Kernel and "Failure" Kernel
        % for the Rank Score
        NBC_KDE_Rank(i,j)  = trapz(x_s, abs(F_s-F_f));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Plot the Kernel Distributions (Less overlap = Smaller Score)
        %figure(i)
        %plot(x_s,F_s,'LineWidth',2);
        %hold on
        %plot(x_f,F_f,'r--','LineWidth',2);
        %title(['For Parameter ', num2str(i), ' in output ', num2str(j)])
        %hold off
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end


%% [7] Rank Parameters and Export Data
% (a) Rank Parameters by sorting the NBC-KDE matrix
toc
for i = 1:n_out
    [score, rank] = sort(NBC_KDE_Rank(:,i), 'descend');
    fprintf('Order of Parameters by NBC Rank for CCU Eval Output %.4f \n', i)
    rank
    sprintf('The Rank Scores for the Above Order are: ')
    NBC_KDE_Rank(:,i)
end

% (b) Export Results
xlswrite('MC-NBC_ParameterSpace.xls', ParamSpace)
xlswrite('MC-NBC_EvaluatedOutputs.xls', outputs)
xlswrite('MC-NBC_EvalMetricTargets.xls', targets)
xlswrite('MC-NBC_ClassifiedOutputs.xls', classmat)
xlswrite('MC-NBC_SuccessKernelData.xls', S_KernelData)
xlswrite('MC-NBC_FailureKernelData.xls', F_KernelData)
xlswrite('MC-NBC_Parameters_Ranked.xls', NBC_KDE_Rank)