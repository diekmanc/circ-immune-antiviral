function z = goyal_gamma_pi_19fold_4_21_25(t,u)

%% immune model

global beta delta k r phi pie gamma omega q delta_E m 

S = u(1);
I = u(2);
V = u(3);
M1 = u(4);
M2 = u(5);
E = u(6);

A = u(7);
C = u(8);
n = u(9);

scaledC = 0.1+1.8*(C+1.13)/2.26;

gamma_circ = gamma * scaledC;

z(1) = (1/24)*(-beta * V * S);
z(2) = (1/24)*(beta * V * S - delta * (I^k)*I - m * ( (E^r)/(E^r + phi^r)) * I);
z(3) = (1/24)*(pie * I - gamma_circ * V);
z(4) = (1/24)*(omega * I * M1 - q * M1);
z(5) = (1/24)*(q * (M1 - M2));
z(6) = (1/24)*(q * M2 - delta_E * E);

%% circadian model

global mu taux kappa lambda eta
global alpha0 zeta0 p G
global tShift lux duty

ft=heaviside(square(2*(pi/24)*(t-tShift),duty));

zeta = lux*ft;

alpha = alpha0*(zeta/zeta0)^p;
Bhat = G*(1-n)*alpha;
B = Bhat*(1-0.4*C)*(1-0.4*A);

z(7) = (pi/12)*(mu*(A-(4/3)*A^3)-C*(((24/(0.99669*taux))^2+kappa*B)));
z(8) = (pi/12)*(A+B);
z(9) = lambda*(alpha*(1-n)-eta*n);

z = z';

end