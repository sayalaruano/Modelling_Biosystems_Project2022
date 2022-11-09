%% PARAMETERS
param.S=10; %implant size
param.R= 0.001; %debris removal rate
param.k0= 1; %M0 migration rate
param.Mmax= 10; %maximal M0
param.I=5; %implant stiffness
param.dM0= 0.1; %M0 apoptosis rate
param.vmax1= 0.1; %MM vmax rate M0 to M1
param.km1= 2; %MM km M0 to M1
param.dM1= 0.05; %M1 apoptosis rate
param.vmax2= 0.1; %MM vmax rate M0 to M2
param.km2= 1; %MM km M0 to M2
param.dM2= 0.01; %M2 apoptosis rate
param.kDc1= 0.01; %constant c1 production rate by debris
param.kM1= 0.1; %constant c1 production rate by M1
param.dc1= 0.5; %c1 removal rate
param.kM2= 0.1; %constant c2 production rate by M2
param.dc2= 0.5; %c2 removal rate 
param.kf= 0.01; %proliferation rate F
param.Fmax= 2; %maximal F
param.df= 0.2; %F apoptosis rate
param.Diffc1=10; %c1 diffusion rate
param.Diffc2=10; %c2 diffusion rate

%% Hill-type parameter estimation
M1_obs = [0; 0.1; 0.3; 0.5;1; 1.9; 3; 3.4; 5; 7.5; 8; 9; 12; 20];
F_obs = [0; 3.04E-11; 7.45E-09; 9.54E-08; 3.01E-06; 7.55E-05; 0.000738; 0.001362; 0.008706; 0.042005; 0.0499; 0.064312; 0.088364; 0.098986];
hn = [1:0.1:10]; % arbitrary values adjusted based on the output of the figure
hk = [0.1:0.01:1]; % k is a rate so it can only take values between 0 and 1

% normalization step
M1_obs = M1_obs/max(M1_obs);
F_obs = F_obs/max(F_obs); 

for i = 1:length(hn)
    for j = 1 :length(hk)
        hill = hill_eq(hk(j),hn(i),M1_obs);
        RMSE(i, j) = sqrt(mean((F_obs - hill).^2));
    end
end


min_RMSE = min(RMSE,[],'all');
[x,y] = find(RMSE==min_RMSE);
param.hk = hk(y); %hill type parameter k (migration rate constant)
param.hn = hn(x); %hill type parameter n (power)

X = linspace(0,1);
hill_final = hill_eq(param.hk, param.hn,X);

figure;
hold on;
plot(M1_obs, F_obs,'*');
plot(X, hill_final);
xlabel('[M1] (type 1 macrophage)')
ylabel('[F] (fibroblast concentration)')
hname = ['Fitted Hill Curve, k = ', num2str(param.hk), ', n = ',num2str(param.hn)];
legend('Observations',hname,'Location','SouthEast')
hold off;


%% EQUATIONS % don't run this code
% change in D over time
vD = -param.R*M2;

% change in C1 over time
vC1 = param.kDc1*D + param.kM1*M1 - param.dc1*C1;

% change in M0 over time
m1 = (param.vmax1*M0)/(param.km1 + M0); %michaelis-menten equation rate for transformation from M0 to M1
m2 = (param.vmax2*M0)/(param.km2 + M0); %michaelis-menten equation rate for transformation from M0 to M2
vM0 = D * param.k0 * (1 -(M0/param.Mmax)) - C1*m1*M0 - C2*m2*M0 - param.dM0*M0;
% change in M1 over time
vM1 = C1*m1*M0 - param.dM1*M0;
% change in M2 over time
vM2 = C2*m2*M0 - param.dM2*M0;
% change in C2 over time
vC2 = param.kM2*M2 - param.dc2*C2;
% change in F over time
vF = param.kf*F*(1-F/param.Fmax) + (M1 / (param.hk + M))^param.hn - param.df*F;


%% INITIAL CONDITION
Debris: (10:50) %(10 - low), %(50 - high) 
M0: 1
C1: 0.1
C2: 0.1 

%%
function hill=hill_eq(k,n,C) % hill type function
    hill = C.^n./(k^n + C.^n);
end
