clear all; close all

%% immune model parameters and initial conditons

global beta delta k r phi pie gamma omega q delta_E m 
global V_initial

% set parameters for patient S18

beta = 10^(-7.22976);
delta = 3.14709;
k = 0.0783677;
pie = 10^(2.59026);
m = 3.19727;
omega = 10^(-4.55828);
q = 2.39E-05;

r = 10;
phi = 100;
delta_E = 1;
gamma = 15;

% set initial conditions

V_initial=pie/gamma;
S_0 = 1e7; I_0 = 1; V_0 = V_initial; M1_0 = 1; M2_0 = 0; E_0 = 0;

%% run simulation

t0 = 0;
tf = 24*30;
tin = 0:(1/1000):tf;

options = odeset('AbsTol',1e-8,'RelTol',1e-8,'Events',@stopGoyal_6_5_24);

inits = [S_0 I_0 V_0 M1_0 M2_0 E_0];

[t,u,te,ye,ie] = ode15s(@goyal_model_6_5_24,tin,inits,options);

V = u(:,3);


%% make plot

close('all')

set(0,'DefaultAxesFontSize',24)

f1=figure(1);
semilogy(t/24,V,'Color','k','linewidth',5)
xlabel('Days','interpreter','latex')
ylabel('Viral Load','interpreter','latex')
set(gca,'box','off','XTick',0:5:30)
xlim([0 25])
f1.Position=[680*1.5 558*1.5 560*1.5 420*1.5];
