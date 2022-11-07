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

%% EQUATIONS % don't run this code
% change in D over time
vD = -param.R*M2;

% change in C1 over time
vC1 = param.kDc1*D + param.KM1*M1 - param.dc1*C1;

% change in M0 over time
m1 = (param.vmax1*M0)/(param.km1 + M0); %michaelis-menten equation rate for transformation from M0 to M1
m2 = (param.vmax2*M0)/(param.km2 + M0); %michaelis-menten equation rate for transformation from M0 to M2
vM0 = D * param.k0 * (1 -(M0/param.Mmax)) - C1*m1*M0 - C2*m2*M0 - param.dM0*M0;
% change in M1 over time
vM1 = C1*m1*M0 - param.dM1*M0;
% change in M2 over time
vM2 = C2*m2*M0 - param.dM2*M0;
% change in C2 over time

% change in F over time
vF = param.kf*F*(1-param.Fmax/F);


%% INITIAL CONDITION
Debris: (10:50) %(10 - low), %(50 - high) 
M0: 1
C1: 0.1
C2: 0.1 
