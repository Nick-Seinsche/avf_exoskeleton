function u = W_interpolate()
    W = zeros(8, 1);
    k = 0;
    alpha = 0;
    h = 2^(-k);
    h_T = h/2;

    wall_incline = tan(alpha / 360 * 2 * pi);
    t11 = [wall_incline; -wall_incline; 1];
    t11 = [-(-1)^(1)*wall_incline; (-1)^(1)*wall_incline; 1];
    t11_v = rot_mat(atan(t11(2)),0,0) * t11;
    z1 = [h_T/2, -t11_v(1)/t11_v(3) * h_T/2];
    z4 = -z1;
    z3 = [h_T/2, h/2];
    z6 = [-h_T/2, h/2];
    z2 = 0.5 * (z3 + z1);
    z5 = 0.5 * (z4 + z6);

    z = [z1', z2', z3', z4', z5', z6'];

    x = sym('x', 'real');
    y = sym('y', 'real');
    t = (y + h_T/2) / h_T;

    temp_Lhat = (x - z1(2)) * (x - z2(2)) * (x - z3(2));
    factorhat = subs(diff(temp_Lhat), z1(2));
    Lhat = temp_Lhat / factorhat;

    temp_Lcheck = (x - z4(2)) * (x - z5(2)) * (x - z6(2));
    factorcheck = subs(diff(temp_Lcheck), z4(2));
    Lcheck = temp_Lcheck / factorcheck;

    L1 = (x - z2(2)) / (z1(2) - z2(2)) * (x - z3(2)) / (z1(2) - z3(2));
    L2 = (x - z1(2)) / (z2(2) - z1(2)) * (x - z3(2)) / (z2(2) - z3(2));
    L3 = (x - z1(2)) / (z3(2) - z1(2)) * (x - z2(2)) / (z3(2) - z2(2));

    L4 = (x - z5(2)) / (z4(2) - z5(2)) * (x - z6(2)) / (z4(2) - z6(2));
    L5 = (x - z4(2)) / (z5(2) - z4(2)) * (x - z6(2)) / (z5(2) - z6(2));
    L6 = (x - z4(2)) / (z6(2) - z4(2)) * (x - z5(2)) / (z6(2) - z5(2));

    factorL1 = subs(diff(L1), z1(2));
    factorL2 = subs(diff(L2), z1(2));
    factorL3 = subs(diff(L3), z1(2));

    factorL4 = subs(diff(L4), z4(2));
    factorL5 = subs(diff(L5), z4(2));
    factorL6 = subs(diff(L6), z4(2));
    
    B1 = L1 - factorL1 * Lhat;
    B2 = L2 - factorL2 * Lhat;
    B3 = L3 - factorL3 * Lhat;

    B4 = L4 - factorL4 * Lcheck;
    B5 = L5 - factorL5 * Lcheck;
    B6 = L6 - factorL6 * Lcheck;

    %plot(linspace(z1(2), z3(2), 100), arrayfun(@(xx) subs(B1, x, xx), linspace(z1(2), z3(2), 100))); hold on
    %plot(linspace(z1(2), z3(2), 100), arrayfun(@(xx) subs(B2, x, xx), linspace(z1(2), z3(2), 100)));
    %plot(linspace(z1(2), z3(2), 100), arrayfun(@(xx) subs(B3, x, xx), linspace(z1(2), z3(2), 100))); 
    %plot(linspace(z1(2), z3(2), 100), arrayfun(@(xx) subs(Lhat, x, xx), linspace(z1(2), z3(2), 100))); hold off
    %set(gca, 'XDir','reverse')
    %legend(["q_1", "q_2", "q_3", "q_4"]);

    max_shift = -t11_v(1)/t11_v(3) * h_T;
    c = max_shift / (h/2 - 2 * (t - 0.5) * max_shift / 2);

    ptop = W(1) * B1 + W(2) * B2 + W(3) * B3 + W(4) * Lhat; 
    pbot = W(5) * B4 + W(6) * B5 + W(7) * B6 + W(8) * Lcheck; 

    g_plus = (1 - t) * (h/2 - x) * c;
    g_minus = -t * (h/2 - x) * c;

    %gg = @(yy,xx) eval(subs(subs((h/2 - x) * c, x, xx), y, yy));

    p = subs(ptop, x, x + g_plus) * t + subs(pbot, x,x + g_minus) * (1-t);

    u = @(yy,xx) eval(subs(subs(p, x, xx), y, yy));
end

function Q = rot_mat(a,b,c)
    Q = [1, 0, 0; 0, cos(a), -sin(a); 0, sin(a), cos(a)];
    Q = Q * [cos(b), 0, sin(b); 0, 1, 0; -sin(b), 0, cos(b)];
    Q = Q * [cos(c), -sin(c), 0; sin(c), cos(c), 0; 0,0,1];
end