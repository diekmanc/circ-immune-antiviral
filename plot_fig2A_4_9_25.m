clear all; close all

%% read in data

% patient S18 is row 11, S12 is row 8, and G2 is 13

inFile1 = sprintf('patientRow%d_circNone_30days_3_7_25.csv',11);
inFile2 = sprintf('patientRow%d_circNone_30days_3_7_25.csv',8);
inFile3 = sprintf('patientRow%d_circNone_30days_3_7_25.csv',13);

data1 = dlmread(inFile1);
data2 = dlmread(inFile2);
data3 = dlmread(inFile3);

t1 = data1(:,1);
t2 = data2(:,1);
t3 = data3(:,1);

V1 = data1(:,4);
V2 = data2(:,4);
V3 = data3(:,4);

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
