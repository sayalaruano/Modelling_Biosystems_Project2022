% Define parameters for the drug treatment
%arbitrary rates
param.drug1 = 0.1; 
param.drug2 = 0.1;
%k for log growth; 
param.kdrug1 = 0.01; 
param.kdrug2 = 0.01; 
% drug1 clearance rate
% assumption: all molecules (like c1 and c2) are removed at the same rate
param.dDrug1 = 0.0005; 
param.dDrug2 = 0.0005;

% Set the initial conditions (nondimensional)
%         D; C1; C2; M0; M1; M2; F; drug1; drug2
x0_low_drug = [10; 0.1; 0.1; 1; 0; 0; 0; 0.2; 0.2]; 
x0_high_drug = [50; 0.1; 0.1; 1; 0; 0; 0; 0.2; 0.2];

% Drug treatment 1 simulation
param.TreatmentOption = 1;
[T_drug1,X_drug1] = ode15s(@fibroblastencaps_drug,[t0 tf],x0_low_drug,[],param);

% Drug treatment 2 simulation
param.TreatmentOption =2;
[T_drug2,X_drug2] = ode15s(@fibroblastencaps_drug,[t0 tf],x0_low_drug,[],param);

% Drug treatment 3 simulation
param.TreatmentOption =3;
[T_bothdrugs,X_bothdrugs] = ode15s(@fibroblastencaps_drug,[t0 tf],x0_low_drug,[],param);

figure
hold on;
% Plot the treatments and base case
plot(T_low,X_low(:,7), 'LineWidth', 1.05);
plot(T_drug1,X_drug1(:,7), 'LineWidth', 1.05);
plot(T_drug2,X_drug2(:,7), 'LineWidth', 1.05);
plot(T_bothdrugs,X_bothdrugs(:,7), 'LineWidth', 1.05);
% Add titles and legend
xlabel('Time');
ylabel('Concentration');
leg = legend('No treatment','C1 inhibition','C2 inhibition','Inhibition of both C1 and C2', 'Location','NorthEast');
title(leg,'Drug treatment'); 
legend('boxoff');
hold off;

function dx = fibroblastencaps_drug(t,x,param)
    D = x(1);
    C1 = x(2);
    C2 = x(3);
    M0 = x(4);
    M1 = x(5);
    M2 = x(6);
    F = x(7);
    drug1 = x(8);
    drug2 = x(9);
    
    % setting drug effect parameters based on treatment option
    param.Kdrug1 = 0; % effect of drug 1 on C1
    param.Kdrug2 = 0; % effect of drug 2 on C2
    if param.TreatmentOption == 1
        param.Kdrug1 = param.drug1; 
    elseif param.TreatmentOption == 2
        param.Kdrug2 = param.drug2;
    elseif param.TreatmentOption == 3
        param.Kdrug1 = param.drug1;
        param.Kdrug2 = param.drug2;
    end
    
    % either increase removal or decrease influence on michaelis menten rates for M1 and M2
    
    ddrug1 =  - param.dDrug1*drug1; %logistic - removal rate
    ddrug2 =  - param.dDrug2*drug2; %logistic - removal rate
    
    % Define equations for the simulation
    % Debris cell pupulation with law of mass action
    dD = -param.R*M2*D; 
    % Michaelis-menten equation rate for transformation from M0 to M1
    m1 = mich_menten(param.vmax1, M0, param.km1);
    dM1 = C1*m1 - param.dM1*M1;
    % Michaelis-menten equation rate for transformation from M0 to M2
    m2 = mich_menten(param.vmax2, M0, param.km2);
    dM2 = C2*m2 - param.dM2*M2;
    % Logistic equation for Macrophages M0
    dM0 = logistic_eq(D, param.Mmax, param.k0, M0) - C1*m1 - C2*m2 - param.dM0*M0; 
    % Pro-inflammatory cytokines (C1) with law of mass action and modified with the drug treatment 
    dC1 = param.kDc1*D + param.kM1*M1 - param.dc1*C1 - param.Kdrug1*drug2*C1;
    % Anti-inflammatory cytokines (C2) with law of mass action and modified with the drug treatment
    dC2 = param.kM2*M2 - param.dc2*C2 - param.Kdrug2*drug2*C2;
    % Fibroblast with the hill equation
    h = hill_eq(param.hk, param.hn, M1, param.FMM);
    dF = logistic_eq(param.kf, param.Fmax, F, F) + h*M1 - param.df*F;
    
    % Return the equations
    dx = [dD; dC1; dC2; dM0; dM1; dM2; dF; ddrug1; ddrug2];
end