
clc
moment = readtable('data_moments.txt'); moment = moment{:,:};
mean_wealth = moment(4);
mean_pi = moment(2);
mean_participation = moment(1);
sR = readmatrix('realized_R.txt');
rt = sR - 1; 
rf = 0.01;
T = 45;  
cap_tax = 0.3;

% Initialize the result variable
result = 0;

for t = 0:T-1
    % Calculate the expressions inside the summation
    expression1 = mean_wealth * (mean_pi * rt(t+1) * (rt(t+1) > 0) + (1 - mean_pi) * rf) * cap_tax / (1 + rf)^t;
    expression2 = mean_wealth * rf * cap_tax / (1 + rf)^t;

    % Update the result variable with the summation
    result = result + mean_participation * expression1 + (1 - mean_participation) * expression2;
end

% Display the final result
disp(result);

