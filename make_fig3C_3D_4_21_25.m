clear all; close all

global beta delta k r phi pie gamma omega q delta_E m 
global EC50
global circFlag treatFlag
global V_initial
global dose
global ka

circFlag = 1;
treatFlag = 0;
dose = 0;
dose1 = 200/70; 
dose2 = 100/70;

params = dlmread('Estimated_params.csv',',',1,1);

% patient S18 is row 11, S5 is row 4, S12 is row 8, S14 is row 9, G2 is 13
ix = 13;

beta = 10^params(ix,2);
delta = params(ix,3);
k = params(ix,4);
pie = 10^params(ix,5);
m = params(ix,6);
omega = 10^params(ix,7);
q = params(ix,10);

r=10;
phi=100;
delta_E = 1;
gamma = 15;
EC50 = 0.8;
ka = 12;

V_initial=pie/gamma;
S_0=1e7;I_0=1; V_0=V_initial; M1_0=1; M2_0=0; E_0=0; Cp_0=0; Ca_0=0;

global mu taux kappa lambda eta
global alpha0 zeta0 p G
global tShift lux duty

photo=12;
lux=1000;
lambda=60;
duty=100*(photo/24);
mu=0.23;
taux=24.2;
kappa=0.55;
alpha0=0.05; eta=0.0075; G=33.75; p=0.5; % 2013 paper
zeta0=9500;

numDays=30;

t0=0;
tf1=24*numDays;
tin1 = t0:(1/1000):tf1;

options = odeset('AbsTol',1e-8,'RelTol',1e-8,'Events',@stopGoyal_4_9_25);

LL=0;
DD=0;

tShift=0;

A_0 = 1.081178209000000;
C_0 = -0.179561130400000;
n_0 = 0.003088693984000;

inits1 = [S_0 I_0 V_0 M1_0 M2_0 E_0 Cp_0 Ca_0 A_0 C_0 n_0];

%% no treatment

t = [];
u = [];
stopFlag = 0;

treatStartDay = 50;

for dx = 1:numDays
    for hx = 1:24
        if dx >= treatStartDay
            treatFlag = 1;
            if hx == treatHour
                if dx == treatStartDay
                    dose = dose1;
                else
                    dose = dose2;
                end
            else 
                dose = 0;
            end
        else 
            treatFlag = 0;
        end
        if stopFlag == 0
            t0 = (hx - 1) + (dx - 1)*24;
            tf = t0 + 1;
            tin = t0:(1/1000):tf;
            if dx == 1 && hx == 1
                inits = inits1;
            else
                inits = u_temp(end,:);
                inits(7) = inits(7) + dose;
            end
            [t_temp,u_temp,te,ye,ie] = ode15s(@goyal_with_k3d_treatment_19fold_4_21_25,tin,inits,options);
            t = [t; t_temp];
            u = [u; u_temp];
            if length(te) > 0
                stopFlag = 1;
            end
        end
    end
end

t1 = t;
I1 = u(:,2);

%% treat at ZT 0 starting on Day 5

t = [];
u = [];
stopFlag = 0;

treatStartDay = 6;
treatHour = 1;

for dx = 1:numDays
    for hx = 1:24
        if dx >= treatStartDay
            treatFlag = 1;
            if hx == treatHour
                if dx == treatStartDay
                    dose = dose1;
                else
                    dose = dose2;
                end
            else 
                dose = 0;
            end
        else 
            treatFlag = 0;
        end
        if stopFlag == 0
            t0 = (hx - 1) + (dx - 1)*24;
            tf = t0 + 1;
            tin = t0:(1/1000):tf;
            if dx == 1 && hx == 1
                inits = inits1;
            else
                inits = u_temp(end,:);
                inits(7) = inits(7) + dose;
            end
            [t_temp,u_temp,te,ye,ie] = ode15s(@goyal_with_k3d_treatment_19fold_4_21_25,tin,inits,options);
            t = [t; t_temp];
            u = [u; u_temp];
            if length(te) > 0
                stopFlag = 1;
            end
        end
    end
end

t2 = t;
I2 = u(:,2);

%% treat at ZT 12 starting on Day 5

t = [];
u = [];
stopFlag = 0;

treatStartDay = 6;
treatHour = 13;

for dx = 1:numDays
    for hx = 1:24
        if dx >= treatStartDay
            treatFlag = 1;
            if hx == treatHour
                if dx == treatStartDay
                    dose = dose1;
                else
                    dose = dose2;
                end
            else 
                dose = 0;
            end
        else 
            treatFlag = 0;
        end
        if stopFlag == 0
            t0 = (hx - 1) + (dx - 1)*24;
            tf = t0 + 1;
            tin = t0:(1/1000):tf;
            if dx == 1 && hx == 1
                inits = inits1;
            else
                inits = u_temp(end,:);
                inits(7) = inits(7) + dose;
            end
            [t_temp,u_temp,te,ye,ie] = ode15s(@goyal_with_k3d_treatment_19fold_4_21_25,tin,inits,options);
            t = [t; t_temp];
            u = [u; u_temp];
            if length(te) > 0
                stopFlag = 1;
            end
        end
    end
end

t3 = t;
I3 = u(:,2);

%% treat at ZT 0 starting on Day 6

t = [];
u = [];
stopFlag = 0;

treatStartDay = 7;
treatHour = 1;

for dx = 1:numDays
    for hx = 1:24
        if dx >= treatStartDay
            treatFlag = 1;
            if hx == treatHour
                if dx == treatStartDay
                    dose = dose1;
                else
                    dose = dose2;
                end
            else 
                dose = 0;
            end
        else 
            treatFlag = 0;
        end
        if stopFlag == 0
            t0 = (hx - 1) + (dx - 1)*24;
            tf = t0 + 1;
            tin = t0:(1/1000):tf;
            if dx == 1 && hx == 1
                inits = inits1;
            else
                inits = u_temp(end,:);
                inits(7) = inits(7) + dose;
            end
            [t_temp,u_temp,te,ye,ie] = ode15s(@goyal_with_k3d_treatment_19fold_4_21_25,tin,inits,options);
            t = [t; t_temp];
            u = [u; u_temp];
            if length(te) > 0
                stopFlag = 1;
            end
        end
    end
end

t4 = t;
I4 = u(:,2);

%% treat at ZT 12 starting on Day 6

t = [];
u=[];
stopFlag = 0;

treatStartDay = 7;
treatHour = 13;

for dx = 1:numDays
    for hx = 1:24
        if dx >= treatStartDay
            treatFlag = 1;
            if hx == treatHour
                if dx == treatStartDay
                    dose = dose1;
                else
                    dose = dose2;
                end
            else 
                dose = 0;
            end
        else 
            treatFlag = 0;
        end
        if stopFlag == 0
            t0 = (hx - 1) + (dx - 1)*24;
            tf = t0 + 1;
            tin = t0:(1/1000):tf;
            if dx == 1 && hx == 1
                inits = inits1;
            else
                inits = u_temp(end,:);
                inits(7) = inits(7) + dose;
            end
            [t_temp,u_temp,te,ye,ie] = ode15s(@goyal_with_k3d_treatment_19fold_4_21_25,tin,inits,options);
            t = [t; t_temp];
            u = [u; u_temp];
            if length(te) > 0
                stopFlag = 1;
            end
        end
    end
end

t5 = t;
I5 = u(:,2);

%% make Figure 3A

set(0,'DefaultAxesFontSize',36)

f1=figure(1);
semilogy(t1/24,I1,'Color','k','linewidth',5)
hold on
semilogy(t2/24,I2,'Color','b','linewidth',5)
semilogy(t3/24,I3,'Color','g','linewidth',5)
semilogy(t4/24,I4,'Color','r','linewidth',5)
semilogy(t5/24,I5,'Color','m','linewidth',5)
semilogy(t1/24,I1,'Color','k','linewidth',5)

xlabel('Days','interpreter','latex')
ylabel('Infected Cells','interpreter','latex')
set(gca,'box','off','XTick',0:5:30)
xlim([0 15])
f1.Position=[680*1.5 558*1.5 560*1.5 420*1.5];
legend('no treatment','Day 5 ZT 0','Day 5 ZT 12', 'Day 6 ZT 0','Day 6 ZT 12','location','northeast','interpreter','latex')
legend('boxoff')
title('G2')


%% make Figure 3B

ind2 = min(find(t2==24*5));
ind3 = min(find(t3==24*5+12));
ind4 = min(find(t4==24*6));
ind5 = min(find(t5==24*6+12));

inds2 = ind2+24024:24024:length(t2);
inds3 = ind3+24024:24024:length(t3);
inds4 = ind4+24024:24024:length(t4);
inds5 = ind5+24024:24024:length(t5);

f2=figure(2);
semilogy(t1/24,I1,'Color','k','linewidth',5)
hold on
semilogy(t2/24,I2,'Color','b','linewidth',5)
semilogy(t3/24,I3,'Color','g','linewidth',5)
semilogy(t4/24,I4,'Color','r','linewidth',5)
semilogy(t5/24,I5,'Color','m','linewidth',5)
semilogy(t1/24,I1,'Color','k','linewidth',5)
semilogy(t2(ind2)/24,I2(ind2),'bo','MarkerSize',20,'MarkerFaceColor','b')
semilogy(t3(ind3)/24,I3(ind3),'go','MarkerSize',20,'MarkerFaceColor','g')
semilogy(t4(ind4)/24,I4(ind4),'ro','MarkerSize',20,'MarkerFaceColor','r')
semilogy(t5(ind5)/24,I5(ind5),'mo','MarkerSize',20,'MarkerFaceColor','m')
semilogy(t2(inds2)/24,I2(inds2),'bo','MarkerSize',20,'linewidth',5)
semilogy(t3(inds3)/24,I3(inds3),'go','MarkerSize',20,'linewidth',5)
semilogy(t4(inds4)/24,I4(inds4),'ro','MarkerSize',20,'linewidth',5)
semilogy(t5(inds5)/24,I5(inds5),'mo','MarkerSize',20,'linewidth',5)

xlabel('Days','interpreter','latex')
ylabel('Infected Cells','interpreter','latex')
set(gca,'box','off','XTick',4:.5:9)
xlim([4 8])
f2.Position=[680*1.5 558*1.5 560*1.5 420*1.5];
title('G2')