% Set the time of the simulation
t0 = 0; tf = 2000;  %(s)

% Set the initial conditions (nondimensional)
%         D; C1; C2; M0; M1; M2; F
x0_low = [10; 0.1; 0.1; 1; 0; 0; 0]; 
x0_high = [50; 0.1; 0.1; 1; 0; 0; 0]; 

% Low initial concentration of debris simulation
[T_low,X_low] = ode15s(@fibroblastencaps,[t0 tf],x0_low,[],param);

% Plot all the low initial concentrations together
figure;
hold on
% Plot all the variables
for i=1:length(x0_low)
    plot(T_low,X_low(:,i), 'LineWidth', 1.05);
end
% Add titles and legend
xlabel('Time');
ylabel('Concentration');
leg = legend('Debris', 'Cytokine 1', 'Cytokine 2','Macrophage 0','Macrophage 1','Macrophage 2','Fibroblast','Location','NorthEast');    
title(leg,'Cell population');
legend('boxoff');
hold off;

% Plot low initial individual concentrations
% Create vectors of names and colors for the cell populations 
names = {['Debris'], ['Cytokine 1'], ['Cytokine 2'], ['Macrophage 0'], ['Macrophage 1'], ['Macrophage 2'], ['Fibroblast'], 'Drug1', 'Drug2'};
options = {[0 0.4470 0.7410], [1 0 0], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840], [1 1 1], [1 1 1]};

% Plot individual variables
for i=1:length(x0_low)
    if i <= 7
        figure
        plot(T_low,X_low(:,i),'Color', options{i}, 'LineWidth', 1.05);
        xlabel('Time');
        ylabel('Concentration');
        leg = legend(names(i),'Location','NorthEast');
        title(leg,'Cell population');
        legend('boxoff');
    end
end

function dx=fibroblastencaps(t,x,param)
    % Assign parameters for the simulation 
    D = x(1);
    C1 = x(2);
    C2 = x(3);
    M0 = x(4);
    M1 = x(5);
    M2 = x(6);
    F = x(7);

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
    % Pro-inflammatory cytokines (C1) with law of mass action
    dC1 = param.kDc1*D + param.kM1*M1 - param.dc1*C1;
    % Anti-inflammatory cytokines (C2) with law of mass action
    dC2 = param.kM2*M2 - param.dc2*C2;
    % Fibroblast with the hill equation
    h = hill_eq(param.hk, param.hn, M1, param.FMM);
    dF = logistic_eq(param.kf, param.Fmax, F, F) + h*M1 - param.df*F;
    
    % Return the equations
    dx = [dD; dC1; dC2; dM0; dM1; dM2; dF];
end

function c =hill_eq(k,n,C,MR) % hill type function
% INPUTS: 
    % k - rate constant
    % n - power 
    % C - independent concentration
    % MR - maximum rate
    c = MR*C.^n./(k^n + C.^n);
end

function m = mich_menten(vmax, concent, km)
% INPUTS: 
    % vmax - maximum rate
    % concent - substrate concentration
    % km - MM constant
    m = (vmax*concent)/(km + concent);
end

function l = logistic_eq(r, K, pop1, pop2)
% INPUTS: 
    % r - growth rate
    % K - carrying capacity 
    % pop1 - population 1
    % pop1 - population 2
    l = r*pop1*(1 -(pop2/K)); 
end