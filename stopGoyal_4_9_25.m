function [value,isterminal,direction] = stopGoyal_4_9_25(t,u)

    value(1) = u(2) - 1;
    
    isterminal(1) = 1;
    
    direction(1) = -1;

end