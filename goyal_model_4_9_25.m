function z = goyal_model_4_9_25(t,u)

global beta delta k r phi pie gamma omega q delta_E m 

S = u(1);
I = u(2);
V = u(3);
M1 = u(4);
M2 = u(5);
E = u(6);

z(1) = (1/24)*(-beta * V * S);
z(2) = (1/24)*(beta * V * S - delta * (I^k)*I - m * ( (E^r)/(E^r + phi^r)) * I);
z(3) = (1/24)*(pie * I - gamma * V);
z(4) = (1/24)*(omega * I * M1 - q * M1);
z(5) = (1/24)*(q * (M1 - M2));
z(6) = (1/24)*(q * M2 - delta_E * E);

z = z';

end