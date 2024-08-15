function u = G_interpolate(V)%, h, h_T, alpha)
    h = 1/2;
    h_T = 1/4;
    alpha = 20.9;
    wall_incline = tan(alpha / 360 * 2 * pi);
    t11 = [wall_incline; -wall_incline; 1];

    z1 = [h_T/2, t11(2)/t11(3) * h_T/2];
    z4 = -z1;
    z3 = [h_T/2, h/2];
    z6 = [-h_T/2, h/2];
    z2 = 0.5 * (z3 + z1);
    z5 = 0.5 * (z4 + z6);

    z = [z1', z2', z3', z4', z5', z6'];
  
    L1 = @(y,x) (y + h_T / 2) / h_T * (x - z2(2)) / (z1(2) - z2(2)) * (x - z3(2)) / (z1(2) - z3(2));
    L2 = @(y,x) (y + h_T / 2) / h_T * (x - z1(2)) / (z2(2) - z1(2)) * (x - z3(2)) / (z2(2) - z3(2));
    L3 = @(y,x) (y + h_T / 2) / h_T * (x - z1(2)) / (z3(2) - z1(2)) * (x - z2(2)) / (z3(2) - z2(2));

    L4 = @(y,x) (1 - (y + h_T / 2) / h_T) * (x - z5(2)) / (z4(2) - z5(2)) * (x - z6(2)) / (z4(2) - z6(2));
    L5 = @(y,x) (1 - (y + h_T / 2) / h_T) * (x - z4(2)) / (z5(2) - z4(2)) * (x - z6(2)) / (z5(2) - z6(2));
    L6 = @(y,x) (1 - (y + h_T / 2) / h_T) * (x - z4(2)) / (z6(2) - z4(2)) * (x - z5(2)) / (z6(2) - z5(2));

    u = @(y,x) V(1) * L1(y,x) + V(2) * L2(y,x) + V(3) * L3(y,x) + V(4) * L4(y,x) + V(5) * L5(y,x) + V(6) * L6(y,x);
    
end