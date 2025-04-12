clear all; close all

%% immune model

global beta delta k r phi pie gamma omega q delta_E m 
global V_initial

% set parameters and initial conditions common to all patients

r = 10;
phi = 100;
delta_E = 1;
gamma = 15;

S_0 = 1e7; I_0 = 1; M1_0 = 1; M2_0 = 0; E_0 = 0; 

% set integration options

options = odeset('AbsTol',1e-8,'RelTol',1e-8,'Events',@stopGoyal_4_9_25);

t0=0;
tf=24*30;
tin = 0:(1/1000):tf;

%% set patient-specific parameters and simulate with no circadian regulation

% patient S18 is row 11, S5 is row 4, S12 is row 8, S14 is row 9, G2 is 13

params = dlmread('Estimated_params.csv',',',1,1);

% patient S18
ix = 11; 

beta = 10^params(ix,2); 
delta = params(ix,3);
k = params(ix,4);
pie = 10^params(ix,5);
m = params(ix,6);
omega = 10^params(ix,7);
q = params(ix,10);

V_initial=pie/gamma;
V_0 = V_initial;
inits = [S_0 I_0 V_0 M1_0 M2_0 E_0];

[t1,u] = ode15s(@goyal_model_6_5_24,tin,inits,options);
V1 = u(:,3);

% patient S12
ix = 8;

beta = 10^params(ix,2); 
delta = params(ix,3);
k = params(ix,4);
pie = 10^params(ix,5);
m = params(ix,6);
omega = 10^params(ix,7);
q = params(ix,10);

V_initial=pie/gamma;
V_0 = V_initial;
inits = [S_0 I_0 V_0 M1_0 M2_0 E_0];

[t2,u] = ode15s(@goyal_model_6_5_24,tin,inits,options);
V2 = u(:,3);

% patient G2
ix = 13; 

beta = 10^params(ix,2); 
delta = params(ix,3);
k = params(ix,4);
pie = 10^params(ix,5);
m = params(ix,6);
omega = 10^params(ix,7);
q = params(ix,10);

V_initial=pie/gamma;
V_0 = V_initial;
inits = [S_0 I_0 V_0 M1_0 M2_0 E_0];

[t3,u] = ode15s(@goyal_model_6_5_24,tin,inits,options);
V3 = u(:,3);


%% make figure

set(0,'DefaultAxesFontSize',24)

figure(1)
semilogy(t1/24,V1,'k','linewidth',5)
hold on
semilogy(t2/24,V2,'r','linewidth',5)
semilogy(t3/24,V3,'g','linewidth',5)
semilogy(t1/24,V1,'k','linewidth',5)
legend('S18','S12','G2')
xlabel('Days','interpreter','latex')
ylabel('Viral Load','interpreter','latex')
set(gca,'box','off','XTick',0:5:30)
xlim([0 25])
legend('boxoff')
