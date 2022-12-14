function vSystem = Fibrosis(InitialStatesVec, param) 

D = InitialStatesVec(1);
M0 = InitialStatesVec(2);
C1 = InitialStatesVec(3);
C2 = InitialStatesVec(4);
M1 = InitialStatesVec(5);
M2 = InitialStatesVec(6);
F = InitialStatesVec(7);

vD = -param.R*M2*D;
m1 = (param.vmax1*M0)/(param.km1 + M0); %michaelis-menten equation rate for transformation from M0 to M1
m2 = (param.vmax2*M0)/(param.km2 + M0); %michaelis-menten equation rate for transformation from M0 to M2
vM0 = D * param.k0 * (1 -(M0/param.Mmax)) - C1*m1 - C2*m2 - param.dM0*M0;
vC1 = param.kDc1*D + param.kM1*M1 - param.dc1*C1;
vC2 = param.kM2*M2 - param.dc2*C2;
vM1 = C1*m1 - param.dM1*M0;
vM2 = C2*m2 - param.dM2*M0;
vF = param.kf*F*(1-F/param.Fmax) + param.fmmax*(M1 / (param.hk + M1))^param.hn - param.df*F;

vSystem = [vD; vM0; vC1; vC2; vM1; vM2; vF];

end
