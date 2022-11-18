% Define measurements of M1 and fibroblast migration rates
M1_obs = [0; 0.1; 0.3; 0.5;1; 1.9; 3; 3.4; 5; 7.5; 8; 9; 12; 20];
F_obs = [0; 3.04E-11; 7.45E-09; 9.54E-08; 3.01E-06; 7.55E-05; 0.000738; 0.001362; 0.008706; 0.042005; 0.0499; 0.064312; 0.088364; 0.098986];

% Set the maximum fibroblast migration rate
param.FMM = max(F_obs);

% Set the vectors of the parameters for calculating the RMSE
hn = 1:0.1:10;
hk = 0.1:0.01:max(M1_obs); 

% Calculate the RMSE values
for i = 1:length(hn)
    for j = 1 :length(hk)
        hill = hill_eq(hk(j),hn(i),M1_obs,param.FMM);
        RMSE(i, j) = sqrt(mean((F_obs - hill).^2));
    end
end

% Find the minimum RMSE value
min_RMSE = min(RMSE,[],'all');
[x,y] = find(RMSE==min_RMSE);

% Find the k (migration rate constant) and n (power) parameters of the minimum RMSE
param.hk = hk(y);
param.hn = hn(x);

% Calculate the Hill type kinetics with the best fitting values
X = linspace(0,max(M1_obs));
hill_final = hill_eq(param.hk, param.hn,X,param.FMM);

% Plot the fitted kinetics and the experimental data
figure;
hold on;
plot(M1_obs, F_obs,'*');
plot(X, hill_final);
xlabel('[M1] (Type 1 macrophage concentration)')
ylabel('[F] (Fibroblast concentration)')
hname = ['Fitted Hill Curve, k = ', num2str(param.hk), ', n = ',num2str(param.hn)];
legend('Observations',hname,'Location','Northwest');
hold off;

function c =hill_eq(k,n,C,MR) % hill type function
% INPUTS: 
    % k - rate constant
    % n - power 
    % C - independent concentration
    % MR - maximum rate
    c = MR*C.^n./(k^n + C.^n);
end