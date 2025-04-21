function z = goyal_with_k3d_treatment_19fold_4_21_25(t,u)

global circFlag treatFlag

global beta delta k r phi pie gamma omega q delta_E m 

global EC50

global ka

S = u(1);
I = u(2);
V = u(3);
M1 = u(4);
M2 = u(5);
E = u(6);
Cp = u(7);
Ca = u(8);

A = u(9);
C = u(10);
n = u(11);

scaledC = 0.1+1.8*(C+1.13)/2.26;

pie_circ = pie * scaledC;

delta_circ = delta;
gamma_circ = gamma;
omega_circ = omega;
beta_circ = beta;
q_circ = q;
delta_E_circ = delta_E;
m_circ = m;
r_circ = r;
phi_circ = phi;
k_circ = k;

kpa = 21;
kc = 29;

V1 = 2.84;
V2 = 0.12;

epsilon = treatFlag * (Ca / V2) / (Ca / V2 + EC50);

z(1) = (1/24)*(-beta_circ * V * S);
z(2) = (1/24)*(beta_circ * V * S - delta_circ * (I^k_circ)*I - m_circ * ( (E^r_circ)/(E^r_circ + phi_circ^r_circ)) * I);
z(3) = (1/24)*(pie_circ * (1 - epsilon) * I - gamma_circ * V);
z(4) = (1/24)*(omega_circ * I * M1 - q_circ * M1);
z(5) = (1/24)*(q_circ * (M1 - M2));
z(6) = (1/24)*(q_circ * M2 - delta_E_circ * E);
z(7) = (1/24)*(-kpa * Cp - kc * Cp);
z(8) = (1/24)*(kpa * Cp - ka * Ca);

%% circadian model

global mu taux kappa lambda eta
global alpha0 zeta0 p G
global tShift lux duty

ft = heaviside(square(2*(pi/24)*(t-tShift),duty));

zeta = lux*ft;

alpha = alpha0*(zeta/zeta0)^p;
Bhat = G*(1-n)*alpha;
B = Bhat*(1-0.4*C)*(1-0.4*A);

z(9) = (pi/12)*(mu*(A-(4/3)*A^3)-C*(((24/(0.99669*taux))^2+kappa*B)));
z(10) = (pi/12)*(A+B);
z(11) = lambda*(alpha*(1-n)-eta*n);

z = z';

end