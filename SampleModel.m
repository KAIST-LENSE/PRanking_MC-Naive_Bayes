function output = SampleModel(kp)             %kp = key parameters
% this is just a trial function for understanding how Sobol Indices work
% with Monte Carlo
x1 = kp(1) ;                 % 1st key parameter
x2 = kp(2) ;                 % 2nd key parameter
x3 = kp(3) ;                 % 3rd key parameter

output = [x1+x2, (x1+2*x3)/x2];      % In this case all parameters are key params
end