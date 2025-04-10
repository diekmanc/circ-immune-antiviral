clear all; close all

%% immune model

global beta delta k r phi pie gamma omega q delta_E m 
global V_initial

% set parameters common to all patients

r = 10;
phi = 100;
delta_E = 1;
gamma = 15;

% set patient-specific parameters
% patient S18 is row 11, S5 is row 4, S12 is row 8, S14 is row 9, G2 is 13

ix = 13;
params = dlmread('Estimated_params.csv',',',1,1);

beta = 10^params(ix,2);
delta = params(ix,3);
k = params(ix,4);
pie = 10^params(ix,5);
m = params(ix,6);
omega = 10^params(ix,7);
q = params(ix,10);

% set initial conditions

V_initial=pie/gamma;
S_0 = 1e7; I_0 = 1; V_0 = V_initial; M1_0 = 1; M2_0 = 0; E_0 = 0; 

inits = [S_0 I_0 V_0 M1_0 M2_0 E_0];

% simulate with no circadian regulation

options = odeset('AbsTol',1e-8,'RelTol',1e-8,'Events',@stopGoyal_3_16_23);

t0=0;
tf=24*30;
tin = 0:(1/1000):tf;

[t0,u] = ode15s(@goyal_model_6_5_24,tin,inits,options);
V0 = u(:,3);

%% circadian model

% set parameters 

global mu taux K
global alpha0 I0 Beta G pp
global tShift lux
global duty
global Phi

photo = 12;
lux = 1000;
Phi = 60;
duty = 100*(photo/24);
mu = 0.23;
taux = 24.2;
K = 0.55;
alpha0 = 0.05; Beta = 0.0075; G = 33.75; pp = 0.5;
I0 = 9500;
tShift=0;

% set initial conditions

A_0 = 1.081178209000000;
C_0 = -0.179561130400000;
n_0 = 0.003088693984000;

inits = [S_0 I_0 V_0 M1_0 M2_0 E_0 A_0 C_0 n_0];

[t1,u] = ode15s(@goyal_circ_beta_19fold,tin,inits,options);
V1 = u(:,3);

[t2,u] = ode15s(@goyal_circ_delta_19fold,tin,inits,options);
V2 = u(:,3);

[t3,u] = ode15s(@goyal_circ_k_19fold,tin,inits,options);
V3 = u(:,3);

[t4,u] = ode15s(@goyal_circ_m_19fold,tin,inits,options);
V4 = u(:,3);

[t5,u] = ode15s(@goyal_circ_r_19fold,tin,inits,options);
V5 = u(:,3);

[t6,u] = ode15s(@goyal_circ_phi,tin,inits,options);
V6 = u(:,3);

[t7,u] = ode15s(@goyal_circ_pi,tin,inits,options);
V7 = u(:,3);

[t8,u] = ode15s(@goyal_circ_gamma,tin,inits,options);
V8 = u(:,3);

[t9,u] = ode15s(@goyal_circ_omega,tin,inits,options);
V9 = u(:,3);

[t10,u] = ode15s(@goyal_circ_q,tin,inits,options);
V10 = u(:,3);

[t11,u] = ode15s(@goyal_circ_delta_E,tin,inits,options);
V11 = u(:,3);

%% make plots

set(0,'DefaultAxesFontSize',36)

C = linspecer(11);

% Figure 2E (all parameters)

f1B=figure(1);
semilogy(t1/24,V1,'Color',C(1,:),'linewidth',5)
hold on
semilogy(t2/24,V2,'Color',C(2,:),'linewidth',5)
semilogy(t3/24,V3,'Color',C(3,:),'linewidth',5)
semilogy(t4/24,V4,'Color',C(4,:),'linewidth',5)
semilogy(t5/24,V5,'Color',C(5,:),'linewidth',5)
semilogy(t6/24,V6,'Color',C(6,:),'linewidth',5)
semilogy(t7/24,V7,'Color',C(7,:),'linewidth',5)
semilogy(t8/24,V8,'Color',C(8,:),'linewidth',5)
semilogy(t9/24,V9,'Color',C(9,:),'linewidth',5)
semilogy(t10/24,V10,'Color',C(10,:),'linewidth',5)
semilogy(t11/24,V11,'Color',C(11,:),'linewidth',5)

legend('$\beta$','$\delta$','$k$','$m$','$r$','$\phi$','$\pi$','$\gamma$','$\omega$','$q$','$\delta_E$','location','EastOutside','interpreter','latex')
legend('boxoff')
f1B.Position=[680*1.5 558*1.5 560*1.5 420*1.5];
xlabel('Days','interpreter','latex')
ylabel('Viral Load','interpreter','latex')
set(gca,'box','off','XTick',0:5:30)
xlim([0 25])
title('G2')